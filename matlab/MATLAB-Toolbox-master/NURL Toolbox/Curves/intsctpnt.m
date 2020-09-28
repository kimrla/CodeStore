function [u, d, pp]=intsctpnt(crv1, crv2 ,u)

% Solve a intersection point of two curves using Newton-Raphson's method
%
%  Inputs: 
% 
%    crv1, crv2 - the two curves to be used
% 
%    u - an initial parametric point of the two curves
%
%  Output: 
% 
%    u - the intersection parametric point
% 
%    d - the minimum distance of the two curves
% 
%    pp - the coordinates of the intersection point of the two curves
%
%  Examples:  
% 
%     p=2; pnt=[6;5.4;0]; 
%     crv1=nrlcirc(3, pnt', pi, 2*pi);
%     x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
%     y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
%     crv2=nrlspline(p,x,y); 
%     uu = aintsctpnts(crv1, crv2);
%     u=uu(3,:)';
%     [u, d, pp]=intsctpnt(crv1, crv2 ,u);
%

% Solve the intersectiont point using Newton-Raphson's method
F=zeros(2,1); dF=zeros(2);
for i=1:10
    [p1, dp1, hessm1] = nrldeval(crv1, u(1)); 
    [p2, dp2, hessm2] = nrldeval(crv2, u(2)); 
    F(1)=2*sum((p1-p2).*dp1); 
    F(2)=-2*sum((p1-p2).*dp2); 
    dF(1,1)=2*sum(dp1.^2)+sum((p1-p2).*(hessm1));
    dF(1,2)=-2*sum(dp1.*dp2);
    dF(2,1)=dF(1,2);
    dF(2,2)=2*sum(dp2.^2)-sum((p1-p2).*(hessm2));
    u=u-dF\F;
    for j=1:2 
        if u(j)<0 
            u(j)=0; 
        elseif u(j)>1 
            u(j)=1; 
        end
    end
end
d=sqrt(sum((p1-p2).^2)); 
pp=[p1,p2];
        
        

%% Demo
% p=2; pnt=[6;5.4;0]; 
% crv1=nrlcirc(3, pnt', pi, 2*pi);
% x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
% y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
% crv2=nrlspline(p,x,y); 
% 
% nrlplot(crv1, 100);
% hold on;
% nrlplot(crv2, 100);
% 
% uu = aintsctpnts(crv1, crv2);
% u=uu(3,:)';
% [u, d, pp]=intsctpnt(crv1, crv2 ,u);
% 
% plot(pp(1,1), pp(2,1), 'ro');
% plot(pp(1,2), pp(2,2), 'bs');





