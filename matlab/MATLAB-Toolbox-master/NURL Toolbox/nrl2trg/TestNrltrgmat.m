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

% Transform the surface into Lagrange blending surface
srf=nrlkntins(srf, {2, 1});
srf = nrldegelev(srf, [2, 2]);
bsrf=nrl2lbf(srf); 
srf=bsrf.faces;

figure; hold on;
nrlplot(srf, [20, 20]);
nrlaxisplot(srf);
view(2); axis equal;

% Coefficients of blending function on triangles
v1 = srf.knots{1}; v2 = srf.knots{2}; 
num1=length(v1); num2=length(v2); 
NB=num1+2*num2-3; 
TN=num1*(num2-1)+1;

% Displacement functions on triangles
[t1, C1]=GaussLobattoR(m, 0, 1); 
[t2, C2]=GaussLobattoR(n, 0, 1);  
t1=linspace(0,1,m); t2=linspace(0,1,n); 
tt={t1, t2}; 

% Basis on an triangle
[G, Gs, Gt, S, T]=nrltrgmat(bsrf, tt);

[s1, r1]=meshgrid(v2(1:end-1), v1(1:end));
pts=[r1(:)', 0; s1(:)', 1];
pnts=nrleval(srf, pts);
p=pnts*G;
dps=pnts*Gs;
dpt=pnts*Gt;

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

% Test for scattered points
[u, v]=meshgrid(tt{1}, tt{2});
[G, Gs, Gt, S, T]=nrltrgmat(bsrf, [u(:)'; v(:)']);
p=pnts*G;
dps=pnts*Gs;
dpt=pnts*Gt;
figure; hold on;
surf(x,y,z);
quiver(p(1,:), p(2,:), dps(1,:), dps(2,:));
quiver(p(1,:), p(2,:), dpt(1,:), dpt(2,:));
colormap summer;      
axis equal;
axis([-0.1, 1.1, -0.1,0.95]);

% k=2;
% Sk=G(k, :); Sk=stretch(Sk,m,n);
% figure; DrawBasis(v1,v2,Sk,S,T);
% str=num2str(k); st=['Basis S', str]; 
% title(st); 
% 
% Skx=Gs(k, :); Skx=stretch(Skx,m,n);
% figure; DrawBasis(v1,v2,Skx,S,T);
% str=num2str(k); st=['Basis S', str, 'x']; 
% title(st);
% 
% Sky=Gt(k, :); Sky=stretch(Sky,m,n);
% figure; DrawBasis(v1,v2,Sky,S,T);
% str=num2str(k); st=['Basis S', str, 'y']; 
% title(st);




