function crv=nrbspline(pts, order)

% Create a nurbs curve using control points of B-splines
% 
% Calling Sequences:
% 
%     crv=nrbspline(pts, order)
%
%  Inputs :   
%         p -    the order of B-spline basis
%        pts =[x; y; z] 
%             x  :   x-coordinates of control points
%             y  :   y-coordinates of control points
%             z  :   z-coordinates of control points
%
%  Example   :  p=2;
%                     x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
%                     y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
%                     crv=nrbspline([x; y], order);
%

[n, m]=size(pts);
if m<n && (m==2 || m==3)
    pts=pts';
    p=m; m=n; n=p;
end
pnts=zeros(3,m);
pnts(1:n,:)=pts;

n=m-1;
L=zeros(n,1);
for i=1:n
    L(i)=sqrt(sum(pnts(:,i).^2+pnts(:,i+1).^2));
end
q=0; 
for i=order+1:n+1
    for j=i-order:i-1
        q=q+L(j);
    end
end
k=n+order+2;
u=zeros(1,k);
for i=order+1:n+1
    s=0;
    for j=i-order:i-1
        s=s+L(j);
    end
    u(i+1)=u(i)+s/q;
end
u(k-order:k)=1;

crv =nrbmak(pnts, u);


%% Demo
% p=2;
% x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
% y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
% crv=nrbspline([x; y], p);
% figure; nrbctrlplot(crv);







