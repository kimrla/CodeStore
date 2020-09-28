function [pnts, ut, tol]=crvtolerance(crv, k, c)

% Get tolerance of length for a crv and sampling points on it
% 
% Calling Sequences:
% 
%     [pnts, tt, tol]=crvtolerance(srf)
%     [pnts, tt, tol]=crvtolerance(srf, k)
%     [pnts, tt, tol]=crvtolerance(srf, k, c)
% 
% INPUTS:
%
%      crv - a nurls curve
%      k - the minimum number of points on each direction 
%           (default value is 13)
%      c - coefficient of the number of sampling points
%           (default value is 1)
%
% OUTPUT:
% 
%     pnts - points on the curve
%     ut - parametric points
%     tol - tolerance of length
% 

if nargin == 1
    k=13; c=1;
end
if nargin == 2
    c=1;
end

[~, cvt]=curvature(crv, crv.knots);
km2=max(abs(cvt));
lnc=nrlmeasure(crv);
k=max([k, c*fix(km2*lnc)]);
tol=lnc/k; 
ut=linspace(0, 1, k);
pnts=nrleval(crv, ut);


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
% u=[]; v=0.5;
% srf=srf1;
% crv=nrlsrf2crv(srf, u, v); 
% nrlplot(crv, 100); 
% view(3);  
% 
% % First and second derivatives of the curve
% [pnts2, ut, tol2]=crvtolerance(crv, 5, 3);
% 
% % First and second derivatives of the surface
% [pnts1, tt, tol1]=srftolerance(srf2);
% 
% % Get approximate intersection points
% [x, pnts1, pnts2, dt]=asrfcrvintersct(pnts1, tt, tol1, pnts2, ut, tol2);
% plot3(pnts1(1,:), pnts1(2,:), pnts1(3,:), 'ro');
% plot3(pnts2(1,:), pnts2(2,:), pnts2(3,:), 'bo');





