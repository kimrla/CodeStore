clear; clc;

% Creat a triangular nurls patch
crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
af=pi/3; a=cos(af); b=sin(af);
crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
crv3 = nrlcirc(1, [0, 0], 0, af);
srf = nrltrgcoons(crv1, crv2, crv3);

% Transform the surface into  a nurls triangle patch
srf=nrlkntins(srf, {2, 1});
srf = nrldegelev(srf, [2, 2]);
tsrf=nrl2trg(srf); 
srf=tsrf.faces;

figure; hold on;
nrlplot(srf, [20, 20]);
nrlaxisplot(srf);
view(2); axis equal;

% Coefficients of blending function on triangles
v1 = srf.knots{1}; v2 = srf.knots{2}; 
num1=length(v1); num2=length(v2); 
NB=num1+2*num2-3; 
TN=num1*(num2-1)+1;

% Basis on an triangle
[G, Gs, Gt, S, T, ~, s, t]=nrlmattrg(tsrf);
m=length(s); n=length(t);

[s1, r1]=meshgrid(v2(1:end-1), v1(1:end));
pts=[r1(:)', 0; s1(:)', 1];
pnts=nrleval(srf, pts);
p=pnts*G';
dps=pnts*Gs';
dpt=pnts*Gt';

% Test for gridded points
figure; hold on;
x=reshape(p(1,:), m, n);
y=reshape(p(2,:), m, n);
z=reshape(p(3,:), m, n);
surf(x,y,z);
quiver(p(1,:), p(2,:), dps(1,:), dps(2,:));
quiver(p(1,:), p(2,:), dpt(1,:), dpt(2,:));
colormap summer;      
axis equal;
axis([-0.1, 1.1, -0.1,0.95]);



