%% ========================================================================
function [overlap, intSurface] = TriangleIntersection2D(V, U, N, ...
  getIntersection, debug)
% Triangles V(V0,V1,V2) and U(U0,U1,U2) are are coplanar. Do they overlap?
% INPUTS:
% N - array(n,3) of surface normals where V(i,:,:) and U(i,:,:) are on the same plane
% V - array(n,3,3) (nFace x 3 dimensions x 3 vertices) of surface #1 vertices
% U - array(n,3,3) (nFace x 3 dimensions x 3 vertices) of surface #2 vertices
%
% OUTPUT:
%   iMsk - N x 1 intersection boolean mask marking which triangles overlap
%   intSurface - intersection surface


%  * parameters: vertices of triangle 1: V0,V1,V2
%  *             vertices of triangle 2: U0,U1,U2
%  * result    : returns 1 if the triangles intersect, otherwise 0

%% Constants needed for creating a mesh based on 3 to 6 points in a circle
tri_mesh{6}  = [1 2 6; 2 4 6; 2 3 4; 4 5 6];
tri_mesh{5}  = [1 2 3; 1 3 4; 4 5 1];
tri_mesh{4}  = [1 2 3; 1 3 4];
tri_mesh{3}  = 1:3;
vertices = [];
faces    = [];
pairs    = [];  % each row corresponds to pair of faces. match row number with face number
nVert    = 0;

%% use edge-edge intersections
overlap = false(size(N,1),1);
i1Idx = [1 1 1 2 2 2 3 3 3];
i2Idx = [3 3 3 1 1 1 2 2 2];
j1Idx = [1 2 3 1 2 3 1 2 3];
j2Idx = [3 1 2 3 1 2 3 1 2];
for row = 1:size(N,1)
  % When it is necesary to project 3D plane on 2D, dIdx will be the optimal
  % dimensions to use.
  [~, a] = max(abs(N(row,:))); 
  [b, c] = otherDim(a); 
  dIdx = [b, c]; 
  order = [];

  %% test all edges of triangle 1 against the edges of triangle 2
  % triangles overlap if edges cross
  [edgeMat, P] = EdgesIntersect3D(...
    squeeze(V(row,:,i1Idx))',squeeze(V(row,:,i2Idx))', ...
    squeeze(U(row,:,j1Idx))',squeeze(U(row,:,j2Idx))');
  overlap(row) = any(edgeMat);
  if ~getIntersection && overlap(row), continue; end
  
  if ~overlap(row)
    %% project onto an axis-aligned plane, that maximizes the area
    % of the triangles, compute indices: dIdx which correspond to 2 smallest N1
    % components.
    V2d = [V(row,dIdx,1); V(row,dIdx,2); V(row,dIdx,3)]; % each row is a 2D vertex
    U2d = [U(row,dIdx,1); U(row,dIdx,2); U(row,dIdx,3)];
    
    %% test if tri1 is totally contained in tri2 or vice varsa
    if PointInTriangle2D(V2d(1,:), U2d) % tri1 is totally contained in tri2
      overlap(row) = true;
      order = 1:3;
    elseif PointInTriangle2D(U2d(1,:), V2d) % tri2 is totally contained in tri1
      overlap(row) = true;
      order = 4:6;
    end
    if overlap(row) && ~getIntersection, continue; end
    clear V2d U2d
  end
  
  %% Build the intersection surface
  if getIntersection && overlap(row)
    %Assemble all the points which might be needed for desining
    %intersection polygon: Intersection points and points from triangle 1
    %and 2
    points   = [P(edgeMat,:); squeeze(V(row,:,:))'; squeeze(U(row,:,:))'];
    if isempty(order) % when one tri is totally contained in the other tri then order is set
      order = IntersectionPolygon(edgeMat>0, points, dIdx, debug);
      if isempty(order), continue; end
    end
    nPoint   = length(order);    % how many points will be added?
    nFace    = nPoint-2;         % how many faces will be added?
    vertices = [vertices; points(order,:)]; %#ok<*AGROW>
    faces    = [faces; nVert+tri_mesh{nPoint} ];
    pairs    = [pairs; row+zeros(nFace,1)];  % each row corresponds to pair of faces. match row number with face number
    nVert    = nVert + nPoint;
    if debug
      assert(max(faces(:))<=size(vertices,1), 'Bad surface definition')
    end
  end
end % for

%% Prepare outputs
intSurface.vertices = vertices;
intSurface.faces    = faces;
if isempty(faces)
  intSurface.edges = [];
else
  intSurface.edges = [faces(:,1:2); faces(:,2:3); faces(:,[1,3])];
end
end % function