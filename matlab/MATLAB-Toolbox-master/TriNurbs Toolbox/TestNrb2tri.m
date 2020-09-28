clear;
% Radius (r) and mesh seed length (h0)
r=1; h0=0.2*r;

% Create a nurbs sphere
circ=nrbline([0,0], [r, 0]);
srf=nrbrevolve(circ, [0,0,0], [0,0,1], pi/2);
% srf=nrbtestsrf;

% Transform a nurbs surface into triangular representation
tnrb=nrb2tri(srf, h0);

% Plot the nodes on parametric domain
figure; hold on;
plot(tnrb.nodes(:,1), tnrb.nodes(:,2), 'ro');
axis equal;
title('Nodes on parametric domain');

% Plot the surface
figure; hold on;
nrbplot(srf, [35, 35]);
view(2); shading interp;

figure; triplot(tnrb.delaunay, tnrb.nodes(:,1), tnrb.nodes(:,2)); 
title('Parametric mesh');

figure; 
trisurf(tnrb.delaunay, tnrb.points(:,1), tnrb.points(:,2), tnrb.points(:,3));
view(2); axis equal;
title('Geometric grid');





