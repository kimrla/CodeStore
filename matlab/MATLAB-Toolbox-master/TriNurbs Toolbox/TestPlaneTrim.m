clear; clc;
close all;
% The mesh seed length (h0)
h0=0.5;

% Create a plane square and a plane cuve
lin1=nrbline([2,1], [9,1]);
lin2=nrbline([2,6], [9,6]);
srf=nrbruled(lin1, lin2);
crv=nrbtestcrv;

% Transform a nurbs surface into triangular representation
tsrf=nrb2tri(srf, h0);
tcrv=nrb2tri(crv, h0);

% The nearest points from the surface to the curve
tol=max([tsrf.seeds(1), tcrv.seeds(1)]);
[p1, p2, d]=nearpnts(tsrf.points, tcrv.points, tol);

% Get the relations from points to triangles of tri-nurbs
p2t1=tnrbpts2tri(tsrf);

% Get the intersection points of a tri-nurbs surface with a curve
[ed1, ed2, pts1, pts2, pnts1, pnts2]=tnrbinterline(tsrf, tcrv, p2t1, p1, p2);

% Plot results
figure; hold on; 
tnrbplot(tsrf);
tnrbplot(tcrv);
plot(pnts1(:,1), pnts1(:,2), 'ro');
plot(pnts2(:,1), pnts2(:,2), 'r*');
axis equal;

figure; hold on;
triplot(tsrf.delaunay, tsrf.nodes(:,1), tsrf.nodes(:,2)); 
plot(pts1(:,1), pts1(:,2), 'k.', 'MarkerSize', 13); 
title('Parametric mesh of the surface'); 
axis equal;





