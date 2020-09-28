%% ========================================================================
function [iMsk, intSurface] = TriangleIntersection3D_Moller(...
  V1, V2, V3, N1, d1, dv, ...
  U1, U2, U3, N2, d2, du, ...
  getIntersection, debug)
%TriangleIntersection3D tests if 2 triangles defined in 3D intersect.
% This is a secondary test following Tomas Moller algorithm
%
% INPUTS:
%   V1, V2, V3, - Nx3 array of surface 1 triangle vertex coordinates
%   U1, U2, U3, - Nx3 array of surface 2 triangle vertex coordinates
%   N1, d1      - Nx3 array of surface 1 triangle plane equations N1.X-d1=0
%   N2, d2      - Nx3 array of surface 2 triangle plane equations N2.X-d2=0
%   dv          - Nx3 array of distances of surface 1 triangle vertices to surface 2 planes
%   du          - Nx3 array of distances of surface 2 triangle vertices to surface 1 planes
%   getIntersection - do we need to output the intersecting surface?
%      Algorithm is much simpler if we do not.
%   debug       - In the debugging mode much more extra "sanity check" test
%      are performed.
%
% OUTPUT:
%   iMsk - N x 1 intersection boolean mask marking which triangles overlap
%   intSurface - intersection surface
%
% ALGORITHM:
% The input triangles are guaranteed to intersect the line of intersection
% of the two planes. Furthermore, these intersections form intervals on
% this line, and the triangles overlap iff these intervals overlap as well.
% Hence, the last part of  the algorithm computes a parametric equation
% L(t) of the line of intersection of the two planes, finds the intervals
% (i.e. scalar intervals on L(t)) for which the line lies inside each
% triangle and performs a one-dimensional interval overlap test.
if debug
  ok = size(N1,2)==3 && size(N2,2)==3 && size(dv,2)==3 && size(du,2)==3 && ...
    size(V1,2)==3 && size(V2,2)==3 && size(V3,2)==3 && ...
    size(U1,2)==3 && size(U2,2)==3 && size(U3,2)==3;
  assert(ok, 'Incorrect array dimensions');
end

%% create strip down versions of MATLAB cross and dot function
cross_prod = @(a,b) [...
  a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
  a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
  a(:,1).*b(:,2)-a(:,2).*b(:,1)];
dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
normalize = @(V) bsxfun(@rdivide,V, sqrt(sum(V.^2,2)));

%% Find intervals of surface 1 and 2 triangles
% compute the scalar intervals on L(t) for which the line lies inside each
% triangle

% Plane creates two open half-spaces. Find the odd vertex, which:
% 1) if no or two vertices are on the plane than pick the vertex which is
%    by itself in its half-space
% 2) if one vertex is on the plane and the other two occupy the same
%    half-space than pick the vertex on the plane
% 3) if one vertex is on the plane and the other two occupy different
%    half-spaces than pick one of the vertices off the plane
% Find vertex using a look-up table "lut" with key calculated based on
% sign of dv and du arrays
lut = [0;3;3;2;1;3;2;2;1;1;2;3;3;0;3;3;2;1;1;2;2;3;1;2;3;3;0];
n = numel(d1);
rows = (1:n)';

%% order surface 1 triangle vertices
a1 = lut(sign(dv)*[9; 3; 1] + 14); % calculate the key and call the look-up table
[b1, c1] = otherDim(a1);
if debug
  assert(all(a1>0), 'Something Wrong: triangles are coplanar')
end
a1 = sub2ind([n,3],rows,a1); % convert row and column IDs to array indecies
b1 = sub2ind([n,3],rows,b1);
c1 = sub2ind([n,3],rows,c1);

%% order surface 2 triangle vertices
a2 = lut(sign(du)*[9; 3; 1] + 14); % calculate the key and call the look-up table
[b2, c2] = otherDim(a2);
if debug
  assert(all(a2>0), 'Something Wrong: triangles are coplanar')
end
a2 = sub2ind([n,3],rows,a2);
b2 = sub2ind([n,3],rows,b2);
c2 = sub2ind([n,3],rows,c2);

%% compute direction of L the line of intersection of 2 planes
% containing 2 triangles. Line L parametric equation: t*D+O=0
D = cross_prod(N1,N2);    % D must be perpendicular to both N1 and N2
[~, maxDim] = max(abs(D),[],2); % compute and index to the largest component of D
if(getIntersection)
  D = normalize(D);
  O = zeros(n,3);
  d = [d1, d2, zeros(n,1)];
  for r =1:n
    N = [N1(r,:); N2(r,:); 0, 0, 0];
    N(3,maxDim(r)) = 1;
    dd = d(r,:)';
    O(r,:) = (N\dd)'; %Solve systems of linear equations N*D3 = d for D3
  end
  clear N d dd
end

%% projection of triangle(V1,V2,V3) and triangle(U1,U2,U3) onto intersection line
% Vp and Up are Nx3 arrays with columns indicating corners of triangles 1 and 2
if(getIntersection)
  Vp=[dot_prod(V1-O,D), dot_prod(V2-O,D), dot_prod(V3-O,D)];
  Up=[dot_prod(U1-O,D), dot_prod(U2-O,D), dot_prod(U3-O,D)];
else
  % Project on one of the axis (closest to the intersection line) instead.
  % Simplified projection is faster and sufficient if we do not need
  % intersection line
  idx = sub2ind([n,3],rows,maxDim);
  Vp = [V1(idx), V2(idx), V3(idx)];
  Up = [U1(idx), U2(idx), U3(idx)];
end
clear V1 V2 V3 U1 U2 U3

%% Calculate surface 1 and 2 triangle intervals
% t1 and t2 are intersection points of surface 1 with the intersection line
% t*D+O=0, and s1 & s2 are intersection points of surface 2 with the same
% line. Tomas Moller algorithm made this section much more complicated
% trying to avoid divisions. However, I could not detect any speed-up.
% Operations (ADD: 12; MUL:4 ; DIV:4 )
t1 = Vp(a1) - (Vp(b1)-Vp(a1)).*dv(a1)./(dv(b1)-dv(a1));
t2 = Vp(a1) - (Vp(c1)-Vp(a1)).*dv(a1)./(dv(c1)-dv(a1));
s1 = Up(a2) - (Up(b2)-Up(a2)).*du(a2)./(du(b2)-du(a2));
s2 = Up(a2) - (Up(c2)-Up(a2)).*du(a2)./(du(c2)-du(a2));

%% Order the intervals as to t1<t2 and s1<s2
msk = t2<t1; % order t1 and t2 so t1<t2
t = t1(msk); t1(msk)=t2(msk); t2(msk)=t; % swap
msk = s2<s1; % order s1 and s2 so s1<s2
t = s1(msk); s1(msk)=s2(msk); s2(msk)=t; % swap

%% Perform THE final test we were preparying for.
% It test for the overlap of 2 1D intervals s1->s2 and t1->t2
iMsk = (s1<t2 & t1<s2);

%% calculate intersection segments
n = nnz(iMsk);
if(getIntersection && n>0)
  % p1 = D*max(t1,s1) + O;    p2 = D*min(t2,s2) + O
  p1 = bsxfun(@times,D(iMsk,:),max(t1(iMsk),s1(iMsk))) + O(iMsk,:);
  p2 = bsxfun(@times,D(iMsk,:),min(t2(iMsk),s2(iMsk))) + O(iMsk,:);
  intSurface.vertices = [p1; p2];
  intSurface.faces    = [1:n; n+1:2*n; n+1:2*n]';
  intSurface.edges    = intSurface.faces(:,1:2);
else
  intSurface.vertices = [];
  intSurface.faces    = [];
  intSurface.edges    = [];
end % if
end % function