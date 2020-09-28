clear; clc;
close all;
% The mesh seed length (h0) and force facter (k0)
h0=1.5; k0=1;

% Create a plane surface
af=0;
crv1=nrbcirc(8, [0,0], af, pi-af);
crv2=nrbcirc(8, [0,0], pi+af, 2*pi-af);
crv2=nrbreverse(crv2);
% srf=nrbruled(crv1, crv2);
srf=nrbrevolve(crv1, [0,0,0], [1,0,0], 2*pi);
% srf=nrbtestsrf;

% Transform a nurbs surface into triangular representation
tsrf=nrb2tri(srf, h0);

figure; hold on; 
tnrbplot(tsrf);
axis equal; view(3);
title('The surface before opimization.');

% Get the boundary of the tri-nurbs surface
p2t2=tnrbpts2tri(tsrf);
bp=triboundary(tsrf, p2t2);
bp=bp{1};

% Optimize the triangles
dt=0.1; it=10;
pfix=sort(RemDuplicate(bp));
tsrf=tnrboptimize(tsrf, p2t2, pfix, dt, it, h0, k0);

% Plot results
figure; hold on; 
tnrbplot(tsrf);
plot3(tsrf.points(bp,1), tsrf.points(bp,2), tsrf.points(bp,3), 'r')
axis equal; view(3);
title('The surface after opimization.');







