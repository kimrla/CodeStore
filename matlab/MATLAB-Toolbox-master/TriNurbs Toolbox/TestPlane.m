

% The mesh seed length (h0) and force facter (k0)
h0=0.4; k0=2;

% Create a plane square and a plane cuve
crv=nrbcirc(2, [5,5], 0, 2*pi);
lin1=nrbline([0,0], [10,0]);
lin2=nrbline([0,10], [10,10]);
srf=nrbruled(lin1, lin2);

% Transform a nurbs surface into triangular representation
tcrv=nrb2tri(crv, 1.2*h0);
tsrf=nrb2tri(srf, h0);

% The nearest points from the surface to the curve
tol=max([tsrf.seeds(1), tcrv.seeds(1)]);
[p1, p2, d]=nearpnts(tsrf.points, tcrv.points, tol);

% Get the relations from points to triangles of tri-nurbs
p2t2=tnrbpts2tri(tsrf);

% Get the intersection points of a curve with a tri-nurbs surface
pnts1=tcrv.points;
[pts2, pnts2, dist]=tnrbpoint2srf(tsrf, tcrv.points);
ed2=tnrbinterline(tsrf, tcrv, p2t2, p1, p2);

% Re-triangulation of nodes close to the trimming curve
[tsrf, pta]=tnrbtrimretri(tsrf, p2t2, ed2, pts2, pnts2);

% Get the boundary of the tri-nurbs surface
p2t2=tnrbpts2tri(tsrf);
bp=triboundary(tsrf, p2t2);
bp=bp{1};

% Optimize the triangles
dt=0.1; it=10;
pfix=sort(RemDuplicate([pta; bp]));
% tsrf=tnrboptimize(tsrf, p2t2, pfix, dt, it, h0, k0);

% Plot results
figure; hold on; 
tnrbplot(tsrf);
tnrbplot(tcrv);
plot3(tsrf.points(pta,1), tsrf.points(pta,2), tsrf.points(pta,3), 'ro')
plot3(tsrf.points(bp,1), tsrf.points(bp,2), tsrf.points(bp,3), 'r')
axis equal;
title('The surface after re-triangulation.');

figure; hold on;
triplot(tsrf.delaunay, tsrf.nodes(:,1), tsrf.nodes(:,2)); 
plot(tsrf.nodes(pta,1), tsrf.nodes(pta,2), 'ro')
plot(pts2(:,1), pts2(:,2));
axis equal;
title('Parametric domian of the surface after re-triangulation.');

