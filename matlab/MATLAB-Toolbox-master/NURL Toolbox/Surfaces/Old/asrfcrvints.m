function [x, pnts1, pnts2, dt]=asrfcrvints(srf, pnts1, tt, crv, pnts2, ut, tol)

% Get approximate intersection points of a curve with a surface
% 
% Calling Sequences:
% 
%     [x, pnts1, pnts2, dt]=asrfcrvintersct(pnts1, tt, tol1, pnts2, ut, tol2)
% 
% INPUTS:
%
%     srf -  a nurls surface
%     pnts1 - points on the surface (see srftolerance) 
%     tt - parametric points {tu tv} of the surface
%     crv - a nurls curve
%     pnts2 - points on the curve (see crvtolerance) 
%     ut - parametric points of the curve
%     tol - tolerance of length
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     pnts1 - points of intersections on the surface
%     pnts2 - points of intersections on the curve
%     dt - distance of the two nearest points
% 

% Get approximate intersection points
[q, p]=apntsints(pnts1, pnts2, tol);

% Get parametric intersection points
[v1, u1]=meshgrid(tt{2}, tt{1}); u1=u1(:); v1=v1(:); 
u1=u1(q)'; v1=v1(q)'; 
x=[u1; v1; ut(p)];

% Further improving accuracy
for j=1:4
    n=length(x(1,:));
    [pnts1, jac1]=nrldeval(srf, x(1:2, :));
    [pnts2, jac2]=nrldeval(crv, x(3, :));
    for i=1:n
        x(:,i)=planelineints(pnts1(:,i), {jac1{1}(:,i), jac1{2}(:,i)}, pnts2(:,i), jac2(:,i), x(:,i));
    end
    pnts1=nrleval(srf, x(1:2, :));
    pnts2=nrleval(crv, x(3, :));
    [q, p]=apntsints(pnts1, pnts2, tol/5);
    pnts1=pnts1(:,q);
    pnts2=pnts2(:,p);
    u1=x(1,q); v1=x(2,q); 
    x=[u1; v1; x(3,p)];

    % Remove duplicated intersections
    [x, uu]=RemDuplicate(x', 1e-3);
    x=x';    
end
pnts1=pnts1(:,uu);
pnts2=pnts2(:,uu);

% Distances
n=length(x(1,:));
dt=zeros(n,1);
for i=1:n
    dt(i)=norm(pnts1(:,i)-pnts2(:,i));
end

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
% % First and second derivatives of the curve
% [pnts2, ut, tol2]=crvtolerance(crv);
% 
% % First and second derivatives of the surface
% [pnts1, tt, tol1]=srftolerance(srf2);
% 
% % Get approximate intersection points
% tol=max([tol1, tol2]); 
% [x, pnts1, pnts2, dt]=asrfcrvints(srf2, pnts1, tt, crv, pnts2, ut, tol);
% 
% plot3(pnts1(1,:), pnts1(2,:), pnts1(3,:), 'ro');
% plot3(pnts2(1,:), pnts2(2,:), pnts2(3,:), 'bo');






