function [x, pnts1, pnts2, dt]=asrfcrvintersct(pnts1, tt, tol1, pnts2, ut, tol2)

% Get approximate intersection points of a curve with a surface
% 
% Calling Sequences:
% 
%     [x, pnts1, pnts2, dt]=asrfcrvintersct(pnts1, tt, tol1, pnts2, ut, tol2)
% 
% INPUTS:
%
%     pnts1 - points on the surface (see srftolerance) 
%     tt - parametric points {tu tv} of the surface
%     tol1 - tolerance of length of the surface
%     pnts2 - points on the curve (see crvtolerance) 
%     ut - parametric points of the curve
%     tol2 - tolerance of length of the curve
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     pnts1 - points of intersections on the surface
%     pnts2 - points of intersections on the curve
%     dt - distance of the two nearest points
% 

% Get approximate intersection points
tol=1.5*max([tol1, tol2]); 
k=length(ut);
dm=DistanceMatrix(pnts2(:,:)', pnts1(:,:)'); 
pp=dm<tol;
qq=zeros(k, 1);
for i=1 : k 
    qi=find(pp(i, :)); 
    if ~isempty(qi) 
        di=DistanceMatrix(pnts2(:,i)', pnts1(:,qi)'); 
        [~, id]=min(di); 
        qq(i)=qi(id); 
    end
end
id=qq~=0; 
p=find(id); q=qq(id); 

% Further selections
mm=length(p);
dt=zeros(mm,1);
for i=1:mm
    dt(i)=dm(p(i), q(i)); 
end
IJ=[];
if mm>1
    ddt=tfdm(1:mm, dt); 
    t=1; 
    for i=1:mm-1 
        if ddt(i)<0 && ddt(i+1)>0 
            mt=[i, i+1];
            [~, I]=min(dt(mt)); 
            IJ(t,1)=mt(I); 
            t=t+1; 
        elseif dt(i)<tol/5
            IJ(t,1)=i; 
            t=t+1; 
        end
    end
end
if isempty(IJ)
    mt=1:mm;
    [~, I]=min(dt(mt));
    IJ=mt(I); 
end
p=p(IJ); q=q(IJ);

% Get parametric intersection points
[v1, u1]=meshgrid(tt{2}, tt{1}); u1=u1(:); v1=v1(:); 
u1=u1(q)'; v1=v1(q)'; 
x=[u1; v1; ut(p)];
dt=dt(IJ);
pnts1=pnts1(:, q); 
pnts2=pnts2(:, p); 


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
% [x, pnts1, pnts2, dt]=asrfcrvintersct(pnts1, tt, tol1, pnts2, ut, tol2);
% plot3(pnts1(1,:), pnts1(2,:), pnts1(3,:), 'ro');
% plot3(pnts2(1,:), pnts2(2,:), pnts2(3,:), 'bo');






