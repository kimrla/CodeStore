function [p1, p2, d]=nearpnts(pnts1, pnts2, tol)

% nearpnts: Get the nearest points between two sets of points
% 
% Calling Sequences:
% 
%       [p1, p2, d]=nearpnts(pnts1, pnts2)
% 
%       [p1, p2, d]=nearpnts(pnts1, pnts2, tol)
% 
% INPUTS:
%
%       pnts1, pnts2 - Two sets of points of size m*dim and m*dim,
%                  where dim means dimension and m, n are the 
%                  number of points.
%
%       tol  -  Tollerance.
%
% OUTPUT:
% 
%       p1, p2 - Indexes of nearest points between two surfaces.
%
%       d - Distances between points on the two surfaces.
%  
% Discriptions:
%   
%        If tol is not included. Then all the nearest point 
%        of pnts2 to pnts1 is returned.
%

dm = DistanceMatrix(pnts1, pnts2);
[d, id]=min(dm, [], 2);
if nargin==2
    p1=1:size(pnts1, 1);
    p2=id;
elseif nargin==3
    p=d<tol;
    p1=find(p);
    p2=id(p);
    d=d(p);
end

%% demo
% % The mesh seed length (h0)
% h0=1.5;
% 
% % Create a nurbs sphere
% L=[7,10,0; 5,1,10];
% center=[5,5,4];
% circ=nrbcirc(4, center, 0, pi);
% srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% 
% % Transform a nurbs surface into quadrangular representation
% qnrb1=nrb2quad(srf1, h0);
% qnrb2=nrb2quad(srf2, h0);
% 
% % The nearest points of the two surfaces
% tol=max([qnrb1.seeds(1), qnrb2.seeds(1)]);
% [p1, p2]=nearpnts(qnrb1.UniqPt, qnrb2.UniqPt, tol);
% 
% % Plot the surface and related quadrangle and line
% figure; hold on;
% trisurf(qnrb1.quad, qnrb1.points(:,1), qnrb1.points(:,2), qnrb1.points(:,3));
% trisurf(qnrb2.tri, qnrb2.points(:,1), qnrb2.points(:,2), qnrb2.points(:,3));
% plot3(qnrb1.UniqPt(p1,1), qnrb1.UniqPt(p1,2), qnrb1.UniqPt(p1,3), 'r*');
% plot3(qnrb2.UniqPt(p2,1), qnrb2.UniqPt(p2,2), qnrb2.UniqPt(p2,3), 'ro');
% view(3); axis equal;




