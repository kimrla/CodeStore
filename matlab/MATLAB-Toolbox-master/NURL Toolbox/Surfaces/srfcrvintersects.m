function [x, pts1, pts2, d, jac1, jac2]=srfcrvintersects(srf, crv, x, it)

% Get intersections of a curve with a surface by Newton-Raphson method
% 
% Calling Sequences:
% 
%     [x, pts1, pts2, d]=srfcrvintersects(srf, crv, x)
% 
% INPUTS:
%
%      srf - a nurls surface
%      crv - a nurls curve
%      x - initial value of parametric intersection point
%      it - the number of iterations
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     nts1 - points of intersections on the surface
%     pts2 - points of intersections on the curve
%     d - distance of the two points
% 

if nargin==3
    it=8;
end

F=zeros(3,1); dF=zeros(3);
for i=1:it
    [pts1, jac1, hess1]=nrldeval(srf, x(1:2)'); 
    [pts2, jac2, hess2]=nrldeval(crv, x(3));  
    dr=pts1-pts2; 
    F(1)=dot(dr, jac1{1}); 
    F(2)=dot(dr, jac1{2}); 
    F(3)=-dot(dr, jac2); 
    dF(1,1)=dot(dr, hess1{1,1})+dot(jac1{1}, jac1{1}); 
    dF(1,2)=dot(dr, hess1{1,2})+dot(jac1{1}, jac1{2}); 
    dF(1,3)=-dot(jac1{1}, jac2); 
    dF(2,1)=dF(1,2); dF(3,1)=dF(1,3); 
    dF(2,2)=dot(dr, hess1{2,2})+dot(jac1{2}, jac1{2}); 
    dF(2,3)=-dot(jac1{2}, jac2); 
    dF(3,2)=dF(2,3); 
    dF(3,3)=-dot(dr, hess2)+dot(jac2, jac2); 
    x=x-dF\F; 
    pp=x>1; x(pp)=1; 
    pp=x<0; x(pp)=0; 
end
d=norm(pts2-pts1);

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
% u=[]; v=0.53;
% srf=srf1;
% crv=nrlsrf2crv(srf, u, v); 
% nrlplot(crv, 100); 
% view(3);  
% 
% % First and second derivatives of the curve
% [pnts2, ut, tol2]=crvtolerance(crv);
% 
% % First and second derivatives of the surface
% [pnts1, tt, tol1]=srftolerance(srf2);
% 
% % Get approximate intersection points
% [xi, pnts1, pnts2, dt]=asrfcrvintersct(pnts1, tt, tol1, pnts2, ut, tol2);
% 
% % Get intersection points by Newton-Raphson method
% k=1; 
% x=xi(:,k);
% [x, pts1, pts2, d]=srfcrvintersects(srf2, crv, x);
% plot3(pts1(1,:), pts1(2,:), pts1(3,:), 'ro');
% plot3(pts2(1,:), pts2(2,:), pts2(3,:), 'bo');








