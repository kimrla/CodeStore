clear; clc;

R=4;  
s1=0; s=2*pi; t1=0; t=pi; 
center=[5, 5, 4]; 
srf2 = nrlsphere(R, center, s1, s, t1, t); 
srf1 = nrltestsrf; 
figure; hold on; 
nrlplot(srf1, [100, 100]); 
nrlplot(srf2, [100, 100]); 
axis equal; view(3); 
shading interp; 

% Get intersection points of two surfaces
m=22; order=3;
[x, pnts1, pnts2, dt]=srfintersct(srf1, srf2, m);
crv=nulpts2crv(pnts1, order);
nrlplot(crv, 100);

% Get the parametric intersection curves
pts1=[x(1,:); x(2,:); zeros(1, length(dt))];
pts2=[x(3,:); x(4,:); zeros(1, length(dt))];
crv1=nulpts2crv(pts1, order);
crv2=nrlmake(pts2, crv1.knots, [0, 1], order);

% Save results
tsrf1=trimsrfmak(srf1, crv, crv1);
tsrf2=trimsrfmak(srf2, crv, crv2);
save trisrfs tsrf1 tsrf2;

% Plot trimmed surfaces
nuv=[1300, 1300];
nt=40;
figure; hold on;
nrltrmplot(crv1, srf1, crv, nuv, nt);

% Plot the results
[~, n]=size(x);
figure; hold on;
nrlplot(crv1, 100);
plot(x(1,:), x(2,:), 'ro');
% for i=1:n
%     text(x(1,i), x(2,i), num2str(i), 'Color', 'black');
% end

figure; hold on;
nrlplot(crv2, 100);
plot(x(3,:), x(4,:), 'ro');
% for i=1:n
%     text(x(3,i), x(4,i), num2str(i), 'Color', 'black');
% end

figure; hold on; 
nrlplot(srf1, [30, 30]); 
% nrlplot(srf1, [20, 20], 'quiver'); 
pc=nrleval(crv, linspace(0,1,50));
plot3(pc(1,:),pc(2,:),pc(3,:),'LineWidth', 1.5);
% nrlplot(crv, 50);
axis equal; view(3); 
shading interp; 

figure; hold on; 
nrlplot(srf2, [30, 30]);
% nrlplot(srf2, [20, 20], 'quiver'); 
nrlplot(crv, 50);
axis equal; view(3); 
shading interp; 




