function [pts2, pnts2, dist]=tnrbpoint2srf(tsrf, pnts1, it)

% tnrbintersect: Get the intersection points of two tri-nurbs surfaces
% 
% Calling Sequences:
% 
%       [pts2, pnts2, dist]=tnrbpoint2srf(tsrf, pnts1)
% 
%       [pts2, pnts2, dist]=tnrbpoint2srf(tsrf, pnts1, it)
% 
% INPUTS:
% 
%       tsrf - Triangular representation of a nurbs surface.
% 
%       pnts1 - The given points on Cartesian coordinates.
%
%       it -  The number of iterations.
%
% OUTPUT: 
%
%       pts2 - Parametric points of the surface.
%
%       pnts2 - The nearest points of the surface to given points.
%
%       dist  -  Distances between points on the surface and the given points.
%

if nargin==2
    it=3; 
end
% i= 0x5f3759df - ( i >> 1 )
% The approximated nearest points to a surface
[~, p2]=nearpnts(pnts1, tsrf.points);

% The nearest points to a surface
np=length(pnts1);
pnts2=zeros(np,3);
pts2=zeros(np,2);
dist=zeros(np,1);
dsrf=nrbderiv(tsrf.nurbs);
for k=1:np
    [pts2(k,:), pnts2(k,:), dist(k)]=nrbsrfreverse(tsrf.nurbs, dsrf, tsrf.nodes(p2(k),:), pnts1(k,:)', it);
end


%% demo
% % The mesh seed length (h0)
% h0=0.8;
% 
% % Create a nurbs surface and a point
% circ=nrbcirc(4, [5,5,4], 0, pi);
% srf1=nrbrevolve(circ, [5,5,4], [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% pts=nrbeval(srf1, {[0.2, 0.4, 0.8]; [0.2, 0.4, 0.8]});
% pnts1=pts(:,:)';
% 
% % Transform a nurbs surface into triangular representation
% tsrf=nrb2tri(srf1, h0);
% 
% % The nearest points of a surface to given points
% [pts2, pnts2, dist]=tnrbpoint2srf(tsrf, pnts1);
% 
% % Plot the results
% figure; hold on;
% tnrbplot(tsrf);
% plot3(pnts1(:,1), pnts1(:,2), pnts1(:,3), 'ro');
% plot3(pnts2(:,1), pnts2(:,2), pnts2(:,3), 'r*');
% axis equal; view(3);
% title('Geometric grid');
% 
% figure; hold on;
% triplot(tsrf.delaunay, tsrf.nodes(:,1), tsrf.nodes(:,2)); 
% plot(pts2(:,1), pts2(:,2), 'ro');
% title('Parametric mesh of the surface'); 
% axis equal;


