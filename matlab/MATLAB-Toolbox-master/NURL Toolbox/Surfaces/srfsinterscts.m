function  [x, pnts1, pnts2, dt]=srfsinterscts(srf1, srf2)

% Get intersection points of two surfaces
% 
% Calling Sequences:
% 
%     [x, pnts1, pnts2, dt]=srfsintersct(srf1, srf2) 
% 
% INPUTS:
%
%     srf1, srf2 - nurls surfaces 
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     pnts1 - points of intersections on surface 1
%     pnts2 - points of intersections on surface 2
%     dt - distance of the two nearest points
% 

% Get intersection points of two surfaces
[x, pnts1, pnts2, dt]=srfsintersct(srf1, srf2);

% Remove duplicated intersections and sort the intersection points
[x, pnts1, pnts2, dt]=optintersects(x, pnts1, pnts2, dt, 0.2/length(dt));

% Get the parametric intersection curves
order=3;
pts1=[x(1,:); x(2,:); zeros(1, length(dt))];
pts2=[x(3,:); x(4,:); zeros(1, length(dt))];
crv1=nulpts2crv(pts1, order);
crv2=nrlmake(pts2, crv1.knots, [0, 1], order);

% Insert more intersection points
[x, pnts1, pnts2, dt]=insertintersects(srf1, srf2, x, pnts1, pnts2, dt, crv1, crv2, pts1);

% Check start and end points
nn=length(dt);
p=find((x<1/nn) & (x~=0));
[x, pnts1, pnts2, dt]=endintersects(srf1, srf2, x, pnts1, pnts2, dt, p);
p=find((1-x<1/nn) & (x~=1));
[x, pnts1, pnts2, dt]=endintersects(srf1, srf2, x, pnts1, pnts2, dt, p);


%% Demo
% R=4;  
% s1=0; s=2*pi; t1=0; t=pi; 
% center=[5, 5, 4]; 
% srf2 = nrlsphere(R, center, s1, s, t1, t); 
% srf1 = nrltestsrf; 
% figure; hold on; 
% nrlplot(srf1, [100, 100]); 
% nrlplot(srf2, [100, 100]);
% axis equal; view(3); 
% shading interp; 
% 
% % Get intersection points of two surfaces
% [x, pnts1, pnts2, dt]=srfsinterscts(srf1, srf2);
% 
% % Get the parametric intersection curves
% nn=length(dt); order=3;
% pts1=[x(1,:); x(2,:); zeros(1, nn)];
% pts2=[x(3,:); x(4,:); zeros(1, nn)];
% crv1=nulpts2crv(pts1, order);
% crv2=nrlmake(pts2, crv1.knots, [0, 1], order);
% 
% % Get intersection points of two surfaces
% crv=nulpts2crv(pnts1, order);
% nrlplot(crv, 100);
% 
% % Plot parametric points
% figure; hold on; 
% nrlplot(crv1);
% for i=1:nn
%     text(x(1,i), x(2,i), num2str(i), 'Color', 'black');
% end
% plot(x(1, :), x(2, :), 'bo');
% 
% figure; hold on; 
% nrlplot(crv2);
% for i=1:nn
%     text(x(3,i), x(4,i), num2str(i), 'Color', 'black');
% end
% plot(x(3, :), x(4, :), 'bo');





