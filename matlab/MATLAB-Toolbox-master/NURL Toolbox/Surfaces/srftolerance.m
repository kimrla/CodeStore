function [pnts, tt, tol, jac]=srftolerance(srf, k, c)

% Get tolerance of length for a surface and sampling points on it
% 
% Calling Sequences:
% 
%     [pnts, tt, tol]=srftolerance(srf)
%     [pnts, tt, tol]=srftolerance(srf, k)
%     [pnts, tt, tol]=srftolerance(srf, k, c)
% 
% INPUTS:
%
%      srf - a nurls surface
%      k - the minimum number of points on each direction 
%           (default value is [13, 14])
%      c - coefficient of the number of sampling points
%           (default value is 1)
%
% OUTPUT:
% 
%     pnts - points on the surface
%     tt - parametric points {tu tv}
%     tol - tolerance of length
% 

if nargin == 1
    k=[13, 14]; c=1;
end
if nargin == 2
    c=1;
end
if length(k)==1
    k=k*ones(2,1);
end

[~, jac, hess]=nrldeval(srf, srf.knots); 
[H, K]=srfcurvatures(jac, hess);
k1=H+real(sqrt(H.^2-K));
km=max(abs(k1));
crvs=nrlextract(srf);
lnu=max([nrlmeasure(crvs(3)), nrlmeasure(crvs(4))]);
lnv=max([nrlmeasure(crvs(1)), nrlmeasure(crvs(2))]);
m=max([k(1), c*fix(km*lnu)]); 
n=max([k(2), c*fix(km*lnv)]); 
s1=linspace(0, 1, m); t1=linspace(0, 1, n);
tt={s1, t1};
if nargout==4
    [pnts, jac]=nrldeval(srf, tt);
else
    pnts=nrleval(srf, tt);
end
tol=max([lnu, lnv])/min([m, n]); 


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





