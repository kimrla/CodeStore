clear; clc;

% Load trimmed surfaces
load trimsrfs tsrf1 tsrf2; 
tsrf=tsrf1;
srf=tsrf.surface; 
crv=tsrf.curve; 
pcrv=tsrf.paras; 

% First derivatives
m=55; n=56; nt=100;
t1=linspace(0, 1, m);
t2=linspace(0, 1, n);
ut=linspace(0, 1, nt);
tt={t1, t2};
[pnts, jac]=nrldeval(srf, tt);
[pts, ja]=nrldeval(crv, ut);

% PLot surface and trimmed curve
figure; hold on;
nrlplot(srf, [100, 100]);
plot3(pts(1,:), pts(2,:), pts(3,:), 'r', 'LineWidth', 0.5);
colormap summer;  
axis equal; view(3); 
shading interp;

% Plot parametric curves
figure; hold on;
trmparaplot(pcrv, [m,n], nt);

ps=nrleval(pcrv, ut);
plot(ps(1,:), ps(2,:), 'ro');
ncrv=nrlmake(ps, ut, [0, 1], 3);
sdata.crvs=ncrv; 
sdata.isline=1; 
save pcrv1 sdata;

% Plot trimmed surface
load brimsrf srfs;
figure; hold on;
for i=1:numel(srfs)-1
    psrfi=srfs(i); 
    psi=nrleval(psrfi, tt);
    ptsi=nrleval(srf, psi(1:2, :));
    x=reshape(ptsi(1,:), m, n);
    y=reshape(ptsi(2,:), m, n);
    z=reshape(ptsi(3,:), m, n);
    surf(x,y,z);
    plot3(x([1,end],:)', y([1,end],:)', z([1,end],:)', 'k');
    plot3(x(:,[1,end]), y(:,[1,end]), z(:,[1,end]), 'k');
end
i=2;
load crownsrf srfs;
srf=tsrf2.surface; 
psrfi=srfs(i); 
psi=nrleval(psrfi, tt);
ptsi=nrleval(srf, psi(1:2, :));
x=reshape(ptsi(1,:), m, n);
y=reshape(ptsi(2,:), m, n);
z=reshape(ptsi(3,:), m, n);
surf(x,y,z);
% plot3(pts(1,:), pts(2,:), pts(3,:), 'b', 'LineWidth', 1.5);
colormap summer;  
axis equal; view(3); 
shading interp;


% Find regions
% [map, pits]=intsmap(tt, ut, pcrv);
% plot(pits(1,:), pits(2,:), 'ro');







