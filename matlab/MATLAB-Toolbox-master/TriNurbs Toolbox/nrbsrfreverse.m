function [pts2, pnts2, dist]=nrbsrfreverse(srf, dsrf, x, pnts1, it)

% tnrbintersect: Get the intersection points of two tri-nurbs surfaces
% 
% Calling Sequences:
% 
%       [pts2, pnts2, dist]=nrbsrfreverse(srf, dsrf, x, pnts1)
% 
%       [pts2, pnts2, dist]=nrbsrfreverse(srf, dsrf, x, pnts1, it)
% 
% INPUTS:
% 
%       srf - A nurbs surface.
% 
%       dsrf - A data structure that represents the first
% 		          derivatives of a NURBS curve, surface or volume.
%
%       x=[u, v]   -  Initial value of the parametric coordinates.
% 
%       pnts1 - A given point on Cartesian coordinate.
%
%       it -  The number of iterations. The default value is 3.
% 
% OUTPUT: 
% 
%       pts2 - Parametric points of the surface.
% 
%       pnts2 - The nearest points of the surface to given points.
% 
%       dist  -  Distances between points on the surface and the given points.
%

if nargin==4
    it=3;
end

x=x(:);
for i=1:it
    [pnt2, jac2]=nrbdeval(srf, dsrf, x(1:2)); 
    dr=pnt2-pnts1; 
    F(1,1)=dot(dr, jac2{1}); 
    F(2,1)=dot(dr, jac2{2}); 
    dF(1,1)=dot(jac2{1}, jac2{1}); 
    dF(1,2)=dot(jac2{1}, jac2{2}); 
    dF(2,1)=dF(1,2); 
    dF(2,2)=dot(jac2{2}, jac2{2}); 
    x=x-dF\F; 
    pp=x>1; x(pp)=1; 
    pp=x<0; x(pp)=0; 
end
pts2=x;
pnts2=nrbeval(srf, pts2); 
dist=norm(pnts2-pnts1);
    
%% demo
% % Create a nurbs surface and a point
% circ=nrbcirc(4, [5,5,4], 0, pi);
% srf1=nrbrevolve(circ, [5,5,4], [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% pts=nrbeval(srf1, {[0.2, 0.4, 0.8]; [0.2, 0.4, 0.8]});
% pnts1=pts(:,:)';
% srf=srf1;
% 
% % The nearest points of a surface to given points
% s=linspace(0,1,10); t=linspace(0,1,11);
% [t,s]=meshgrid(t,s);
% pts=[s(:)'; t(:)'];
% points=nrbeval(srf, pts)';
% [~, p2]=nearpnts(pnts1, points);
% np=length(pnts1);
% pnts2=zeros(np,3);
% pts2=zeros(np,2);
% dist=zeros(np,1);
% dsrf=nrbderiv(srf);
% for k=1:np
%     [pts2(k,:), pnts2(k,:), dist(k)]=nrbsrfreverse(srf, dsrf, pts(:, p2(k)), pnts1(k,:)');
% end
% 
% % Plot the results
% figure; hold on;
% nrbplot(srf, [20, 20]);
% plot3(pnts1(:,1), pnts1(:,2), pnts1(:,3), 'ro');
% plot3(pnts2(:,1), pnts2(:,2), pnts2(:,3), 'r*');
% axis equal; view(3);
% title('Geometric grid');
    
    
    
    
    
    