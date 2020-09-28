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
%       pnts1, pnts2 - Two sets of points of size m*dim and n*dim,
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
    p=d<1.618*tol;
    p1=find(p);
    p2=id(p);
    d=d(p);
end

%% demo
% % The mesh seed length (h0)
% h0=0.8;
% 
% % Create a nurbs sphere
% circ=nrbcirc(4, [5,5,4], 0, pi);
% srf1=nrbrevolve(circ, [5,5,4], [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% 
% % Transform a nurbs surface into triangular representation
% tnrb1=nrb2tri(srf1, h0);
% tnrb2=nrb2tri(srf2, h0);
% 
% % The nearest points of two surfaces
% tol=max([tnrb1.seeds(1), tnrb2.seeds(1)]);
% [p1, p2, d]=nearpnts(tnrb1.points, tnrb2.points, tol);
% 
% % Plot the results
% figure; hold on;
% tnrbplot(tnrb1);
% tnrbplot(tnrb2);
% axis equal; view(3);
% title('Geometric grid');
% plot3(tnrb1.points(p1,1), tnrb1.points(p1,2), tnrb1.points(p1,3), 'ro');
% plot3(tnrb2.points(p2,1), tnrb2.points(p2,2), tnrb2.points(p2,3), 'r*');




