clear; clc;

% The number integration nodes (m, n)
m=15; n=16;

% Creat a triangular nurls patch
crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
a=sin(pi/4);
crv2 = nrlline([0.0 0.0 0.0]',[a a 0.0]');
crv3 = nrlcirc(1, [0, 0], 0, pi/4);
srf = nrltrgcoons(crv1, crv2, crv3);

% Rearrange its direction and plot
[srf, du]=trgsrfdirect(srf);

figure; hold on;
nrlplot(srf, [20, 20]);
nrlaxisplot(srf);
view(2); axis equal;

% Quiver plot
figure; hold on;
nrlplot(srf, [20, 20], 'ctrl');
nrlplot(srf, [20, 20], 'quiver');







