clear; clc;

% Give a point or a set of points
pts=[1,1,5; 9,1,9; 9,9,11; 1,9,6; 5,5,8]';

% Mesh seed length
h0=2; 

% Create a nurbs surface
% srf=nrltestsrf;
circ=nrlcirc(5, [5,5,0], 0, pi);
srf=nrlrevolve(circ, [5,5,0], [1,0,0]);
figure; nrlplot(srf, [20, 20]);
hold on; 
shading interp;
plot3(pts(1,:), pts(2,:), pts(3,:), 'ro');

% Get the nerest points
[d, uvm, pntm]=anearestpts(srf, pts, h0);

plot3(pntm(1,:), pntm(2,:), pntm(3,:), 'r*');





