%% ========================================================================
function inside = PointInTriangle2D(V1, U)
% check if V1 is inside triangle U (U1,U2,U3)
% Algorithm is checking on which side of the half-plane created by the
% edges the point is. It uses sign of determinant to calculate orientation
% of point triplets.
% INPUTS:
%   V1 - 1 x 2 coordinates of a point
%   U  - 3 x 2 coordinates of endpoints of 3 edges of a triangle
% OUTPUT:
%   inside - a boolean or boolean array
det2 = @(A,B,C) (A(:,1)-C(:,1))*(B(:,2)-C(:,2)) - (B(:,1)-C(:,1))*(A(:,2)-C(:,2));
b1 = (det2(U(1,:), U(2,:), V1) > 0);
b2 = (det2(U(2,:), U(3,:), V1) > 0);
b3 = (det2(U(3,:), U(1,:), V1) > 0);
inside = ((b1 == b2) & (b2 == b3)); % inside if same orientation for all 3 edges
end % function