function [x, pnts1, pnts2, dt]=srfcrvintersct(srf, crv, pnts1, tt, tol1)

% Get intersection points of a curve with a surface
% 
% Calling Sequences:
% 
%     [x, pnts1, pnts2, dt]=srfcrvintersct(srf, crv, pnts1, tt, tol1)
% 
% INPUTS:
%
%     srf - a nurls surface
%     crv - a nurls curve
%     pnts1 - points on the surface (see srftolerance) 
%     tt - parametric points {tu tv} of the surface
%     tol1 - tolerance of length of the surface
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     pnts1 - points of intersections on the surface
%     pnts2 - points of intersections on the curve
%     dt - distance of the two nearest points
% 

% First and second derivatives of the curve
[pnts2, ut, tol2]=crvtolerance(crv);

% Get approximate intersection points
tol=max([tol1, tol2]); 
[x, pnts1, pnts2, dt]=asrfcrvints(srf, pnts1, tt, crv, pnts2, ut, tol);

% Get intersection points by Newton-Raphson method
tol=max([tol1, tol2]);
n=length(x(1,:)); dt=zeros(n,1);
for i=1:n
    [x(:,i), pnts1(:,i), pnts2(:,i), dt(i)]=srfcrvintersects(srf, crv, x(:,i), 1);
end
dc=dt<tol*1e-10;
x=x(:,dc); 
pnts1=pnts1(:,dc); 
pnts2=pnts2(:,dc); 
dt=dt(dc);


%% Demo
% R=4;  
% s1=0; s=2*pi; t1=0; t=pi; 
% center=[5, 5, 4]; 
% srf1 = nrlsphere(R, center, s1, s, t1, t); 
% srf2 = nrltestsrf; 
% figure; hold on; 
% nrlplot(srf1, [100, 100]); 
% nrlplot(srf2, [100, 100]); 
% axis equal; view(3); 
% shading interp; 
% 
% % Extract a curve from a surface
% u=[]; v=0.38;
% srf=srf1;
% crv=nrlsrf2crv(srf, u, v); 
% nrlplot(crv, 100); 
% view(3);  
% 
% % First and second derivatives of the surface
% [pnts1, tt, tol1]=srftolerance(srf2);
% 
% % Get intersection points by Newton-Raphson method
% [x, pnts1, pnts2, dt]=srfcrvintersct(srf2, crv, pnts1, tt, tol1);
% 
% plot3(pnts1(1,:), pnts1(2,:), pnts1(3,:), 'ro');
% plot3(pnts2(1,:), pnts2(2,:), pnts2(3,:), 'bo');




