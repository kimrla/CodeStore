function [pts, x, d]=interpoint2quad(Q, P, x)

% interline2quad: Get nrearest point of a point to a quadrangle
% 
% Calling Sequences:
% 
%     [pts1, x, d, inter]=interline2quad(Q, L, x)
% 
% INPUTS:
%
%      Q=[P1x, P1y, P1z
%           P2x, P2y, P2z
%           P3x, P3y, P3z
%           P4x, P4y, P4z]
%           Three vertexes of a quadrangle.
%
%      P=[Px, Py, Pz], the point.
% 
%      x =[u; v], initial guess parametric points of 
%           intersection for the quadrangle.
%           The default values are [0.5; 0.5];
% 
% OUTPUT:
% 
%     pts - Points of intersections on the quadrangle.
% 
%     x =[u, v], where [u, v] are parametric points of 
%          intersection for the quadrangle.
% 
%     d - Distance of the two points
% 

if nargin==2
    x=[0.5; 0.5];
end
if islinetri(Q)
    error('The triangle is a straight line!');
end
x=x(:);
% Use twice Newton iterative method to calculate the par-cor
for j=1:2
    [pts, jac]=quadpoint(Q, x(1), x(2));
    dr=pts-P(:)'; 
    F(1,1)=dot(dr, jac{1}); 
    F(2,1)=dot(dr, jac{2}); 
    dF(1,1)=dot(jac{1}, jac{1}); 
    dF(1,2)=dot(jac{1}, jac{2}); 
    dF(2,1)=dF(1,2);
    dF(2,2)=dot(jac{2}, jac{2}); 
    x(1:2)=x(1:2)-dF\F; 
    pts=quadpoint(Q, x(1), x(2)); 
    d=norm(P-pts);
end


%% demo
% % The points
% P=[0.7,0.4,0.6];
% Q=[0,0,0; 1,0,1; 1,1,0; 0,1,1];
% 
% % Get the inter section of the line to the quandrangle
% [pq, x]=interpoint2quad(Q, P);
% 
% % Test of quadpoint
% m=5; n=6;
% s=linspace(0,1,m);
% t=linspace(0,1,n);
% [t, s]=meshgrid(t, s);
% [Pt, Jac, Hess]=quadpoint(Q, s, t);
% Px=reshape(Pt(:,1), m, n);
% Py=reshape(Pt(:,2), m, n);
% Pz=reshape(Pt(:,3), m, n);
% 
% % Plot results
% figure; surf(Px, Py, Pz);
% hold on;
% plot3(P(:,1), P(:,2), P(:,3), 'm*');
% plot3(pq(:,1), pq(:,2), pq(:,3), 'ro');
% axis equal;




