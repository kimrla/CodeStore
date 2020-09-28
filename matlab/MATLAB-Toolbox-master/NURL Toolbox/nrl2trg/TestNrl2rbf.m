% Test triangular element
clear; clc;

% The number integration nodes (m, n)
m=10; n=12; 

% Creat a triangular nurls patch
crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
af=pi/3; a=cos(af); b=sin(af);
crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
crv3 = nrlcirc(1, [0, 0], 0, af);
srf = nrltrgcoons(crv1, crv2, crv3);
srf=nrlkntins(srf, {5, 4});
srf = nrldegelev(srf, [2, 2]);

% Transform nurls surface into Lagrange blending surface
bsrf=nrl2lbf(srf);

% Plot the nurls surface
figure; hold on;
nrlaxisplot(bsrf.faces);
view(2); axis equal;

% Evaluate points and first derivatives of the lbf surface
figure;
lbfplot(bsrf, [m, n], 'quiver');
view(2);






    






