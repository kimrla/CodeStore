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
[x, pnts1, pnts2, dt]=srfsinterscts(srf1, srf2);

% Get the parametric intersection curves
nn=length(dt); order=3;
pts1=[x(1,:); x(2,:); zeros(1, nn)];
pts2=[x(3,:); x(4,:); zeros(1, nn)];
crv1=nulpts2crv(pts1, order);
crv2=nrlmake(pts2, crv1.knots, [0, 1], order);

% Get intersection points of two surfaces
crv=nulpts2crv(pnts1, order);
nrlplot(crv, 100);

% Save results
tsrf1=trimsrfmak(srf1, crv, crv1);
tsrf2=trimsrfmak(srf2, crv, crv2);
save trimsrfs tsrf1 tsrf2;

% Plot parametric points
figure; hold on; 
nrlplot(crv1);
% for i=1:nn
%     text(x(1,i), x(2,i), num2str(i), 'Color', 'black');
% end
% plot(x(1, :), x(2, :), 'bo');
box on;
axis equal;
axis([0,1,0,1]);

figure; hold on; 
nrlplot(crv2);
% for i=1:nn
%     text(x(3,i), x(4,i), num2str(i), 'Color', 'black');
% end
% plot(x(3, :), x(4, :), 'bo');
box on;
axis equal;
axis([0,1,0,1]);





