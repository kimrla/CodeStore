%% ========================================================================
function overlap = TriangleIntersection3D_Rapid( v1, v2, v3, u1, u2, u3, n1, n2 )
%TriangleIntersection3D tests if 2 triangles defined in 3D intersect.
%
% INPUTS:
%   v1, v2, v3, - Nx3 array of surface 1 triangle vertex coordinates
%   u1, u2, u3, - Nx3 array of surface 2 triangle vertex coordinates
%   n1, n2      - Nx3 array of surface 1 & 2 triangle plane normals. Those
%      are optional and if provided than the first 2 steps of the algorithm
%      (which are equivalent to first 2 steps of Moller algorithm) will be 
%      skipped.
%
% OUTPUT:
%   iMsk - N x 1 intersection boolean mask marking which triangles overlap
%
% ALGORITHM:
%   translated from the UNC-CH V-Collide RAPID code
%    https://wwwx.cs.unc.edu/~geom/papers/COLLISION/vcol.pdf

global V1 V2 V3 U1 U2 U3

cross_prod = @(a,b) [...
  a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
  a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
  a(:,1).*b(:,2)-a(:,2).*b(:,1)];

%% shift t1 and t2 by p1#
V1 = zeros(size(v1));
V2 = v2-v1;
V3 = v3-v1;
U1 = u1-v1;
U2 = u2-v1;
U3 = u3-v1;
clear v1 v2 v3 u1 u2 u3

if(nargin<7)
  %% now begin the series of tests
  n1 = cross_prod( V2-V1, V3-V1 ); % face normals
  n2 = cross_prod( U2-U1, U3-U1 ); % face normals
end
  
%% test the face normals
overlap = project6(n1) & project6(n2);
V1 = V1(overlap,:);
V2 = V2(overlap,:);
V3 = V3(overlap,:);
U1 = U1(overlap,:);
U2 = U2(overlap,:);
U3 = U3(overlap,:);
n1 = n1(overlap,:);
n2 = n2(overlap,:);

%% compute triangle edges
e1 = V2-V1;
e2 = V3-V2;
e3 = V1-V3;
f1 = U2-U1;
f2 = U3-U2;
f3 = U1-U3;

%% run more tests
overlap2 = project6(cross_prod(e1, f1));
overlap2 = project6(cross_prod(e1, f2)) & overlap2;
overlap2 = project6(cross_prod(e1, f3)) & overlap2;
overlap2 = project6(cross_prod(e2, f1)) & overlap2;
overlap2 = project6(cross_prod(e2, f2)) & overlap2;
overlap2 = project6(cross_prod(e2, f3)) & overlap2;
overlap2 = project6(cross_prod(e3, f1)) & overlap2;
overlap2 = project6(cross_prod(e3, f2)) & overlap2;
overlap2 = project6(cross_prod(e3, f3)) & overlap2;
overlap2 = project6(cross_prod(e1, n1)) & overlap2;
overlap2 = project6(cross_prod(e2, n1)) & overlap2;
overlap2 = project6(cross_prod(e3, n1)) & overlap2;
overlap2 = project6(cross_prod(f1, n2)) & overlap2;
overlap2 = project6(cross_prod(f2, n2)) & overlap2;
overlap2 = project6(cross_prod(f3, n2)) & overlap2;
overlap(overlap) = overlap2;
end