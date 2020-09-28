%% ========================================================================
function [intersect, X] = EdgesIntersect3D(V1,V2, U1,U2)
%EdgesIntersectPoint3D calculates point of intersection of 2 coplanar
% segments in 3D
%
% INPUTS:
%   V1,V2 - 1 x 3 coordinates of endpoints of edge 1
%   U1,U2 - 1 x 3 coordinates of endpoints of edge 2
% OUTPUT:
%   X - 1 x 3 coordinates of the intersection point
A = V2-V1;
B = U1-U2;
C = U1-V1;
%% Solve system of equations [A,B,1] * [d;e;0] = C for d and e
det3 = @(a,b) ... % determinant of a matrix with columns: [a, b, 1]
  a(:,1).*b(:,2)-a(:,3).*b(:,2) + ...
  a(:,2).*b(:,3)-a(:,2).*b(:,1) + ...
  a(:,3).*b(:,1)-a(:,1).*b(:,3);
f=det3(A,B); % https://en.wikipedia.org/wiki/Cramer%27s_rule#Explicit_formulas_for_small_systems
t=det3(C,B)./f; % use Cramer's rule
s=det3(A,C)./f;
intersect = (t>=0 & t<=1 & s>=0 & s<=1);
X = V1 + bsxfun(@times,A,t);
end % function