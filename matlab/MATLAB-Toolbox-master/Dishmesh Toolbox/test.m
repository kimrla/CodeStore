h0=0.2;

% Create a nurbs sphere
circ=nrlcirc(1, [0,0], 0, pi);
srf=nrlrevolve(circ, [0,0,0], [1,0,0], 2*pi);
% srf=nrltestsrf;

% Transform a nurbs surface into triangular representation
tnrl=nrl2tri(srf, h0);

% Plot the nodes on parametric domain
figure; hold on;
plot(tnrl.nodes(:,1), tnrl.nodes(:,2), 'ro');
axis equal;
title('Nodes on parametric domain');

% Plot the surface
figure; hold on;
nrlplot(srf, [35, 35]);
view(3);
shading interp;

figure; triplot(tnrl.delaunay, tnrl.nodes(:,1), tnrl.nodes(:,2)); 
title('Parametric mesh');

figure; 
trisurf(tnrl.delaunay, tnrl.points(:,1), tnrl.points(:,2), tnrl.points(:,3));
axis equal;
title('Geometric grid');