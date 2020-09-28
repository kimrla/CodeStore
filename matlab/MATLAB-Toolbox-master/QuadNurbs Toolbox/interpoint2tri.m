function [pts, x, d]=interpoint2tri(T, P, x)

% interpoint2tri: Get the nearest point on a triangle to a point P
% 
% Calling Sequences:
% 
%     [pts, x, d]=interpoint2tri(T, P, x)
% 
% INPUTS:
%
%      T=[P1x, P1y, P1z
%           P2x, P2y, P2z
%           P3x, P3y, P3z]
%           Three vertexes of a triangle.
%
%      P=[Px, Py, Pz]
%           A point.
% 
%      x =[u, v], initial guess parametric points of 
%           nearest point on the triangle to a point P.
%           The default values are [0.5; 0.5];
% 
% OUTPUT:
% 
%     pts - Nearest point on a triangle to the point P.
% 
%     x =[u, v], parametric points of nearest point on 
%          the triangle to the point P.
% 
%     d - Distance of the point to the triangle.
% 

if nargin==2
    x=[0.5; 0.5];
end
t=islinetri(T);
if t
    error('The triangle is a straight line!');
end
x=x(:);
[pts, jac]=tripoint(T(1,:), T(2,:), T(3,:), x(1), x(2));
dr=pts-P; 
F(1,1)=dot(dr, jac{1}); 
F(2,1)=dot(dr, jac{2}); 
dF(1,1)=dot(jac{1}, jac{1}); 
dF(1,2)=dot(jac{1}, jac{2}); 
dF(2,1)=dF(1,2);
dF(2,2)=dot(jac{2}, jac{2}); 
x=x-dF\F; 
pts=tripoint(T(1,:), T(2,:), T(3,:), x(1), x(2));
if nargout>2
    d=norm(P-pts);
end

%% demo
% % Create a triangle and a point
% T=[0,0,0; 1,0,0; 1,1,1];
% P=[0.5,0.5,0];
% % P=tripoint(T(1,:), T(2,:), T(3,:), 0.3, 0.3);
% 
% % Get the nrearest point on a triangle to a point
% [pts, x, d]=interpoint2tri(T, P);
% 
% % Plot the results
% n=10;
% tri=tridelaunay(n);
% [u, v]=TrigNodeVect(n);
% Ps=tripoint(T(1,:), T(2,:), T(3,:), u, v);
% figure; hold on;
%  trisurf(tri, Ps(:,1),Ps(:,2),Ps(:,3));
% plot3(P(1), P(2), P(3), 'ro');
% plot3(pts(1), pts(2), pts(3), 'r*');
% axis equal; view(3); 




