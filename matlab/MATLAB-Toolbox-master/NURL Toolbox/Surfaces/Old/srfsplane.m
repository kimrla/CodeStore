function [x, pnts]=srfsplane(srf1, srf2, x, tol)

% Get intersection of two surfaces by tangent planes
% 
% Calling Sequences:
%
%     [x, pnts]=srfsplane(srf1, srf2, x, tol)
%
% INPUTS:
%
%      srf1  - surfaces 1
%      srf2  - surfaces 2
%      x - initial guess of the intersections
%      tol - tolerace
%
% OUTPUT:
%
%   x - parametric coordinates of the intersections
%   pnts - coordinates of the intersections
%

% Remove duplicated intersections and sort the intersection points
dt=zeros(length(x(1,:)),1); 
[x, ~, ~, dt]=optintersects(x, dt', dt', dt);
n=length(dt);
dp=zeros(1, n-1); 
for i=1:n-1
    dp(i)=norm(x(1:2, i+1)-x(1:2, i));
end
ap=sum(dp)/(n-1);
p=dp<2*ap;
p=[true, p];
x=x(:,p); 

% Approximate by plane intersections
[pnts1, jac1]=nrldeval(srf1, x(1:2,:)); 
[pnts2, jac2]=nrldeval(srf2, x(3:4,:)); 
[E1, F1, G1, g1, nor1]=quadratic1(jac1);
[E2, F2, G2, g2, nor2]=quadratic1(jac2);
D10=sum(nor1.*pnts1);
D20=sum(nor2.*pnts2);
[~, n]=size(pnts1);
qq=1:n;
pnts=pnts1;
for i=1:3
    [pnts(:,qq), qq]=srf2plane(pnts1(:,qq), nor1(:,qq), nor2(:,qq), D10(:,qq), D20(:,qq), tol, i);    
    if isempty(qq)
        break;
    end
end
dp1=pnts-pnts1; dp2=pnts-pnts2;
dE1=dot(dp1, jac1{1}(:,:));
dG1=dot(dp1, jac1{2}(:,:));
dE2=dot(dp2, jac2{1}(:,:));
dG2=dot(dp2, jac2{2}(:,:));
for i=1:n
    x(1:2,i)=x(1:2,i)+[G1(i), -F1(i); -F1(i), E1(i)]*[dE1(i); dG1(i)]/g1(i);
    x(3:4,i)=x(3:4,i)+[G2(i), -F2(i); -F2(i), E2(i)]*[dE2(i); dG2(i)]/g2(i);
end

[x, pnts, ~, ~]=optintersects(x, pnts, dt', dt);
pp=x>1 | x<0;
n=length(x(1,:));
q=true(1,n);
for i=1:n
    if ~isempty(find(pp(:,i), 1))
        q(i)=false;
    end
end
x=x(:,q);
pnts=pnts(:,q);


%% Demo
% R=4;  
% s1=0; s2=2*pi; t1=0; t2=pi; 
% center=[5, 5, 4]; 
% srf1 = nrlsphere(R, center, s1, s2, t1, t2); 
% srf2 = nrltestsrf; 
% % srf2=nrl4surf([10.0 0.0 5.3], [0.0 0.0 3.3], [0.0 10.0 5.3], [10.0 10.0 3.3]);
% figure; hold on; 
% nrlplot(srf1, [100, 100]); 
% nrlplot(srf2, [100, 100]); 
% axis equal; view(3); 
% shading interp; 
% 
% % Tolerance of length and sampling points on surface 1
% [pnts1, tt1, tol1]=srftolerance(srf1, [31, 32]);
% [pnts2, tt2, tol2]=srftolerance(srf2, [31, 32]);
% tol=max([tol1, tol2]);
% 
% % Get approximate intersection points
% dm=DistanceMatrix(pnts1(:,:)', pnts2(:,:)'); 
% mn=size(pnts1); mn=mn(2)*mn(3);
% pp=dm<tol;
% qq=zeros(mn, 1);
% for i=1:mn
%     qi=find(pp(i,:));
%     if ~isempty(qi)
%         di=DistanceMatrix(pnts1(:,i)', pnts2(:,qi)'); 
%         [~, id]=min(di);
%         qq(i)=qi(id);
%     end
% end
% id=qq~=0;
% p=find(id); q=qq(id);
% 
% % Get parametric points by plane intersections
% [v1, u1]=meshgrid(tt1{2}, tt1{1}); u1=u1(:); v1=v1(:); 
% [v2, u2]=meshgrid(tt2{2}, tt2{1}); u2=u2(:); v2=v2(:); 
% u1=u1(p); v1=v1(p); 
% u2=u2(q); v2=v2(q); 
% x=[u1'; v1'; u2'; v2']; 
% [x, pnts]=srfsplane(srf1, srf2, x, tol);
% [pnts, uu]=RemDuplicate(pnts', tol/5);
% pnts=pnts';
% x=x(:, uu);
% 
% % Further improving accuracy
% n=length(x(1,:)); dt=zeros(n,1);
% for i=1:n
%     [x(1:2,i), pnts(:,i), dt(i)]=distpntsrf(srf1, pnts(:,i)', x(1,i), x(2,i), 5);
%     [x(3:4,i), pnts(:,i), dt(i)]=distpntsrf(srf2, pnts(:,i)', x(3,i), x(4,i), 5);  
% end
% p=dt<tol/10;
% x=x(:,p); dt=dt(p);
% pnts=pnts(:,p);
% [x, pnts, ~, dt]=optintersects(x, pnts, dt', dt);
% 
% % Get the parametric intersection curves
% order=2;
% n=length(x(1,:));
% pts1=[x(1,:); x(2,:); zeros(1, n)];
% pts2=[x(3,:); x(4,:); zeros(1, n)];
% crv1=nulpts2crv(pts1, order);
% crv2=nrlmake(pts2, crv1.knots, [0, 1], order);
% crv=nulpts2crv(pnts, order);
% 
% plot3(pnts(1,:), pnts(2,:), pnts(3,:), 'ro');
% nrlplot(crv, 100); 
% % Plot results
% figure; hold on; 
% nrlplot(srf1, [100, 100]); 
% plot3(pnts(1,:), pnts(2,:), pnts(3,:), 'ro');
% nrlplot(crv, 100); 
% shading interp; 
% axis equal; view(3); 
% 
% figure; hold on; 
% nrlplot(srf2, [100, 100]); 
% plot3(pnts(1,:), pnts(2,:), pnts(3,:), 'ro');
% nrlplot(crv, 100); 
% shading interp; 
% axis equal; view(3); 
% 
% n=length(x(1,:));
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







