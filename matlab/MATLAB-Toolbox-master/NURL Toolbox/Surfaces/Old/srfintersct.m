function  [x, pnts1, pnts2, dt]=srfintersct(srf1, srf2, m)

% Get intersection points of two surfaces
% 
% Calling Sequences:
% 
%     [x, pnts1, pnts2, dt, pnts, tt, tol]=srfintersct(srf1, srf2, m)
% 
% INPUTS:
%
%     srf1, srf2 - nurls surfaces
%     m - the number of curves to be extracted from srf2
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     pnts1 - points of intersections on surface 1
%     pnts2 - points of intersections on surface 2
%     dt - distance of the two nearest points
% 

% Tolerance of length and sampling points on a surface
[pnts, tt, tol]=srftolerance(srf1);

% Get intersection points of two surfaces
x=[]; pnts1=[]; pnts2=[]; dt=[];
v=linspace(0, 1, m);
for i=1:m
    % Extract a curve from a surface    
    crv=nrlsrf2crv(srf2, [], v(i)); 

    % Get intersection points of the curve with the surface
    [xi, pnts1i, pnts2i, dti]=srfcrvintersct(srf1, crv, pnts, tt, tol);
    xi(4,:)=v(i);
    
    % Assemble the results
    x=[x, xi]; dt=[dt; dti];
    pnts1=[pnts1, pnts1i];
    pnts2=[pnts2, pnts2i];
    
    % Extract a curve from a surface
    crv=nrlsrf2crv(srf2, v(i), []); 

    % Get intersection points of the curve with the surface
    [xi, pnts1i, pnts2i, dti]=srfcrvintersct(srf1, crv, pnts, tt, tol);
    xi(4,:)=xi(3,:); xi(3,:)=v(i);
    
    % Assemble the results
    x=[x, xi]; dt=[dt; dti];
    pnts1=[pnts1, pnts1i];
    pnts2=[pnts2, pnts2i];
end

% Remove duplicated intersections and sort the intersection points
[x, pnts1, pnts2, dt]=optintersects(x, pnts1, pnts2, dt);


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
% m=22; order=3;
% [x, pnts1, pnts2, dt]=srfintersct(srf1, srf2, m);
% crv=nulpts2crv(pnts1, order);
% nrlplot(crv, 100);
% 
% % Get the parametric intersection curves
% pts1=[x(1,:); x(2,:); zeros(1, length(dt))];
% pts2=[x(3,:); x(4,:); zeros(1, length(dt))];
% crv1=nulpts2crv(pts1, order);
% crv2=nrlmake(pts2, crv1.knots, [0, 1], order);
% 
% % Save results
% tsrf1=trimsrfmak(srf1, crv, crv1);
% tsrf2=trimsrfmak(srf2, crv, crv2);
% save trisrfs tsrf1 tsrf2;
% 
% % Plot trimmed surfaces
% nuv=[1300, 1300];
% nt=40;
% figure; hold on;
% nrltrmplot(crv1, srf1, crv, nuv, nt);
% 
% % Plot the results
% [~, n]=size(x);
% figure; hold on;
% nrlplot(crv1, 100);
% plot(x(1,:), x(2,:), 'ro');
% for i=1:n
%     text(x(1,i), x(2,i), num2str(i), 'Color', 'black');
% end
% 
% figure; hold on;
% nrlplot(crv2, 100);
% plot(x(3,:), x(4,:), 'ro');
% for i=1:n
%     text(x(3,i), x(4,i), num2str(i), 'Color', 'black');
% end
% 
% figure; hold on; 
% nrlplot(srf1, [30, 30]); 
% nrlplot(srf1, [20, 20], 'quiver'); 
% nrlplot(crv, 50, 'quiver');
% axis equal; view(3); 
% shading interp; 
% 
% figure; hold on; 
% nrlplot(srf2, [30, 30]);
% nrlplot(srf2, [20, 20], 'quiver'); 
% nrlplot(crv, 50, 'quiver');
% axis equal; view(3); 
% shading interp; 


