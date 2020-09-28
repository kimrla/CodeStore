function [intMatrix, intSurface] = SurfaceIntersection(surface1, surface2, varargin)
%SURFACEINTERSECTION intersection of 2 surfaces
% [intMatrix, intSurface] = SurfaceIntersection(surface1, surface2)
% calculates the intersection of surfaces 1 and 2. Code can either return
% just the matrix indicating which face of surface1 intersected with face
% of surface2, which is calculated using Tomas Moller algorithm, or can
% also return the actual line of intersection. In case when parts of the
% surface 1 and 2 lay on the same plane the intersection is a 2D area
% instead of 1D edge. In such a case the intersection area will be
% triangulated and intSurface.edges will hold the edges of the
% triangulation surface and intSurface.faces will hold the faces.
%
% INPUT:
%  * surface1 & surface2 - two surfaces defined as structs or classes.
%    Several inputs are possible:
%    - struct with "faces" and "vertices" fields
%    - 'triangulation' class (only the boundary surface will be used)
%    - 'delaunayTriangulation' class
%
% OUTPUT:
% * intMatrix - sparse Matrix with n1 x n2 dimension where n1 and n2 are
%               number of faces in surfaces
% * intSurface - a structure with following fields:
%     intSurface.vertices - N x 3 array of unique points
%     intSurface.edges    - N x 2 array of edge vertex ID's
%     intSurface.faces    - N x 3 array of face vertex ID's
%
% ALGORITHM:
% Based on Triangle/triangle intersection test routine by Tomas Mller, 1997.
%  See article "A Fast Triangle-Triangle Intersection Test",
%  Journal of Graphics Tools, 2(2), 1997
%  http://web.stanford.edu/class/cs277/resources/papers/Moller1997b.pdf
%  http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/opttritri.txt

%% Get FACES and VERTICES inputs
if isa(surface1, 'triangulation')
  [surface1.faces, surface1.vertices] = freeBoundary(surface1);
elseif isa(surface1, 'delaunayTriangulation')
  S = surface1;
  surface1 = [];
  surface1.faces    = S.ConnectivityList;
  surface1.vertices = S.Points;
  clear S
end
if isa(surface2, 'triangulation')
  [surface2.faces, surface1.vertices] = freeBoundary(surface2);
elseif isa(surface2, 'delaunayTriangulation')
  S = surface2;
  surface2 = [];
  surface2.faces    = S.ConnectivityList;
  surface2.vertices = S.Points;
  clear S
end
ok1 = isstruct(surface1) && isfield(surface1, 'vertices') && isfield(surface1, 'faces');
ok2 = isstruct(surface2) && isfield(surface2, 'vertices') && isfield(surface2, 'faces');
assert(ok1, 'Surface #1 must be a struct with "faces" and "vertices" fields' );
assert(ok2, 'Surface #2 must be a struct with "faces" and "vertices" fields' );

%% Flip dimentions if necessery
if size(surface1.faces,1)==3 && size(surface1.faces,2)~=3
  surface1.faces = surface1.faces';
end
if size(surface1.vertices,1)==3 && size(surface1.vertices,2)~=3
  surface1.vertices = surface1.vertices';
end
if size(surface2.faces,1)==3 && size(surface2.faces,2)~=3
  surface2.faces = surface2.faces';
end
if size(surface2.vertices,1)==3 && size(surface2.vertices,2)~=3
  surface2.vertices = surface2.vertices';
end

%% Parse extra parameters
getIntersection = (nargout>1);
debug = true;
PointRoundingTol = 1e6;
algorithm = 'moller';
k=1;
nVarargs = length(varargin);
while (k<=nVarargs)
  assert(ischar(varargin{k}), 'Incorrect input parameters')
  switch lower(varargin{k})
    case 'debug'
      debug = varargin{k+1}~=0;
      k = k+1;
    case 'algorithm'
      algorithm = lower(strtrim(varargin{k+1}));
      k = k+1;
    case 'pointroundingtol'
      PointRoundingTol = varargin{k+1};
      k = k+1;
  end
  k = k+1;
end

%% Initialize variables
epsilon = eps;
nFace1 = size(surface1.faces,1);
nFace2 = size(surface2.faces,1);
nVert1 = size(surface1.vertices,1);
nVert2 = size(surface2.vertices,1);

%% create strip down versions of MATLAB cross and dot function
cross_prod = @(a,b) [...
  a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
  a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
  a(:,1).*b(:,2)-a(:,2).*b(:,1)];
dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
normalize = @(V) bsxfun(@rdivide,V, sqrt(sum(V.^2,2)));

%% Initialize output variables
% intersect is a nFace1 x nFace2 matrix. Possible values: -2 (do not know),
% -1 (coplanar with unknown overlap), 0 (no intersections), 1 (intersects).
% Negative values are internal only.
intMatrix  = zeros([nFace1,nFace2], 'int8')-2; % -2 indicates that there was no succesful test yet
intSurface.vertices = [];
intSurface.faces    = [];
intSurface.edges    = [];

% =======================================================================
%% === Stage 1 ==========================================================
% =======================================================================
% Each triangle is a subset of the plane it lies in, so for two triangles
% to intersect they must overlap along the line of intersection of their
% planes. Hence, a necessary condition for intersection is that each
% triangle must intersect the plane of the other.
% Mllers method begins by checking the mutual intersection of each
% triangle with the plane of the other. To do so, it determines for each
% triangle on which side of the other triangles supporting plane its
% vertices lie. Now, if all vertices of one triangle lie on the same side
% and no vertex is on the plane, the intersection is rejected.

%% compute plane equations for each triangle of the surface #1
% plane equation #1: N1.X-d1=0
V1 = surface1.vertices(surface1.faces(:,1),:);
V2 = surface1.vertices(surface1.faces(:,2),:);
V3 = surface1.vertices(surface1.faces(:,3),:);
N1 = cross_prod(V2-V1,V3-V1); % array size nFace1 x 3
N1 = normalize(N1);
d1 = dot_prod(N1,V1);         % array size nFace1 x 1

%% Distance from surface #2 vertices to planes of surface #1
% Calculate signed distance from all vertices of surface #2 to each plane
% of of surface #1
du = zeros(nFace1,nVert2);
for iVert2 = 1:nVert2
  p = surface2.vertices(iVert2,:);
  du(:,iVert2) = N1(:,1)*p(1) + N1(:,2)*p(2) + N1(:,3)*p(3) - d1;
end
if debug
  assert(all(size(du)==[nFace1,nVert2]), 'Incorrect array dimensions: dv')
end
du(abs(du)<epsilon)=0; % robustness check
% Distances from vertex 1, 2 & 3 of faces of surface #2 to planes of surface #1
du1 = du(:,surface2.faces(:,1));
du2 = du(:,surface2.faces(:,2));
du3 = du(:,surface2.faces(:,3));
if debug
  assert(all(size(du1)==size(intMatrix)), 'Incorrect array dimensions: du1')
end
clear du
intMatrix(du1.*du2>0 & du1.*du3>0) = 0;   % same sign on all of them & not equal 0
if(all(intMatrix==0)), return; end        % no intersections
intMatrix(du1==0 & du2==0 & du3==0) = -1; % coplanar with unknown overlap

%% compute plane of triangle (U0,U1,U2)
% plane equation 2: N2.X-d2=0
U1 = surface2.vertices(surface2.faces(:,1),:);
U2 = surface2.vertices(surface2.faces(:,2),:);
U3 = surface2.vertices(surface2.faces(:,3),:);
N2 = cross_prod(U2-U1,U3-U1); % array size nFace1 x 3
N2 = normalize(N2);
d2 = dot_prod(N2,U1);        % array size nFace1 x 1

%% Distance from surface #1 vertices to planes of surface #2
% Calculate signed distance from all vertices of surface #1 to each plane
% of of surface #2
dv = zeros(nFace2,nVert1);
for iVert1 = 1:nVert1
  p = surface1.vertices(iVert1,:);
  dv(:,iVert1) = N2(:,1)*p(1) + N2(:,2)*p(2) + N2(:,3)*p(3) - d2;
end
if debug
  assert(all(size(dv)==[nFace2,nVert1]), 'Incorrect array dimensions: dv')
end
dv(abs(dv)<epsilon)=0; % robustness check
% Distances from vertex 1, 2 & 3 of faces of surface #1 to planes of surface #2
dv1 = dv(:,surface1.faces(:,1))';
dv2 = dv(:,surface1.faces(:,2))';
dv3 = dv(:,surface1.faces(:,3))';
if debug
  assert(all(size(dv1)==size(intMatrix)), 'Incorrect array dimensions: dv1')
end
clear dv
intMatrix(dv1.*dv2>0 & dv1.*dv3>0) = 0;   % same sign on all of them & not equal 0
if(all(intMatrix==0)), return; end        % no intersections
intMatrix(dv1==0 & dv2==0 & dv3==0) = -1; % coplanar with unknown overlap

% =======================================================================
%% === Stage 2 ==========================================================
% =======================================================================

%% Process remaining (non-coplanar) triangle pairs
tMsk = (intMatrix==-2);
n = nnz(tMsk);
if n>0
  [face1, face2] = find(tMsk);
  switch lower(algorithm)
    case 'moller'
      if size(dv1(tMsk),1)==1
        dv = [dv1(tMsk)', dv2(tMsk)', dv3(tMsk)'];
        du = [du1(tMsk)', du2(tMsk)', du3(tMsk)'];
      else
        dv = [dv1(tMsk), dv2(tMsk), dv3(tMsk)];
        du = [du1(tMsk), du2(tMsk), du3(tMsk)];
      end
      
      [intMatrix(tMsk), intSurface] = TriangleIntersection3D_Moller(...
        V1(face1,:), V2(face1,:), V3(face1,:), N1(face1,:), d1(face1,:), dv, ...
        U1(face2,:), U2(face2,:), U3(face2,:), N2(face2,:), d2(face2,:), du, ...
        getIntersection, debug);
    case 'rapid'
      % Undocumented experimental feature. In some experiments I got
      % identical results as with Moller algorithm, but others gave
      % different results. Often faster tham Moller.
      intMatrix(tMsk) = TriangleIntersection3D_Rapid( ...
        V1(face1,:), V2(face1,:), V3(face1,:), ...
        U1(face2,:), U2(face2,:), U3(face2,:), N1(face1,:), N2(face2,:) );
    otherwise
      error('Unknown algorithm name');
  end
end % if

%% Process coplanar triangle pairs. Pass #1:
% compare the overlap of the bounding boxes
tMsk = (intMatrix==-1);
if nnz(tMsk)>0
  [face1, face2] = find(tMsk);
  overlap = true;
  for idim = 1:3
    v = [V1(face1,idim), V2(face1,idim), V3(face1,idim)];
    u = [U1(face2,idim), U2(face2,idim), U3(face2,idim)];
    t1 = min(v,[],2);
    t2 = max(v,[],2);
    s1 = min(u,[],2);
    s2 = max(u,[],2);
    overlap = overlap & (s1<=t2 & t1<=s2);
  end
  % if overlap intMatrix will remain "-1" otherwise it will change to "0"
  intMatrix(tMsk) = -1*overlap;
  clear v u t1 t2 s1 s2 overlap
end

%% Process coplanar triangle pairs. Pass #2:
% use edge-edge intersections
tMsk = (intMatrix==-1);
if nnz(tMsk)>0
  [face1, face2] = find(tMsk);
  
  % repack data prior to function call
  V(:,:,1)=V1(face1,:); V(:,:,2)=V2(face1,:); V(:,:,3)=V3(face1,:);
  U(:,:,1)=U1(face2,:); U(:,:,2)=U2(face2,:); U(:,:,3)=U3(face2,:);
  [intMatrix(tMsk), intSurface2] = TriangleIntersection2D(V, U, ...
    N1(face1,:), getIntersection, debug);
  
  %% Merge surfaces
  if getIntersection
    np = size(intSurface.vertices,1);
    intSurface.vertices = [intSurface.vertices; intSurface2.vertices];
    intSurface.faces    = [intSurface.faces;    intSurface2.faces+np];
    intSurface.edges    = [intSurface.edges;    intSurface2.edges+np];
    if debug
      np = size(intSurface.vertices,1);
      assert(max(intSurface.faces(:))<=np, 'Bad surface definition')
      assert(max(intSurface.edges(:))<=np, 'Bad surface definition')
    end
  end
end

%% Clean up the outputs
intMatrix = sparse(double(intMatrix));
if(getIntersection)
  % make point array unique
  P = round(intSurface.vertices*PointRoundingTol)/PointRoundingTol;
  [~,ia,ic] = unique(P,'rows'); % V = P(ia,:) and P = V(ic,:).
  intSurface.vertices = intSurface.vertices(ia,:);
  intSurface.faces = ic(intSurface.faces);
  intSurface.edges = ic(intSurface.edges);
end
end % function


















