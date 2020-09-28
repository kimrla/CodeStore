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

% Evaluate the derivatives of an triganle with respect to area coodinates
t1=linspace(0,1,m); t2=linspace(0,1,n); 
tt={t1, t2}; 
[pnt, jac]=nrltrgdeval(bsrf, tt);
dps=jac{1}; dpt=jac{2};

figure; hold on;
x=squeeze(pnt(1,:,:));
y=squeeze(pnt(2,:,:));
z=squeeze(pnt(3,:,:));
surf(x,y,z);
quiver(pnt(1,:), pnt(2,:), dps(1,:), dps(2,:));
quiver(pnt(1,:), pnt(2,:), dpt(1,:), dpt(2,:));
colormap summer;      
axis equal;
axis([-0.1, 1.1, -0.1,0.95]);

[u, v]=meshgrid(tt{1}, tt{2});
[pnt, jac]=nrltrgdeval(bsrf, [u(:)'; v(:)']);
dps=jac{1}; dpt=jac{2};
figure; hold on;
surf(x,y,z);
quiver(pnt(1,:), pnt(2,:), dps(1,:), dps(2,:));
quiver(pnt(1,:), pnt(2,:), dpt(1,:), dpt(2,:));
colormap summer;      
axis equal;
axis([-0.1, 1.1, -0.1,0.95]);







