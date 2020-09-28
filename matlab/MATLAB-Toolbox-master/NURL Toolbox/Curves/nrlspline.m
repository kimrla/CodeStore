function crv=nrlspline(p, x, y, z)

% Create a nurl curve using control points of B-splines
% 
% Calling Sequences:
% 
%     crv=nrlspline(p, x, y)
%     crv=nrlspline(p, x, y, z)
%
%  Inputs :   p, x, y, z
%                 p  :    the order of B-spline basis
%                 x  :   x-coordinates of control points
%                 y  :   y-coordinates of control points
%                 z  :   z-coordinates of control points
%
%  Example   :  p=2;
%                     x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
%                     y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
%                     crv=nrlspline(p, x, y);
%

if nargin==3
    z=0*x;
end

n=length(x);
pnts = zeros(3,n);
pnts(1,:)=x(:);
pnts(2,:)=y(:);
pnts(3,:)=z(:);
n=n-1;
L=zeros(n,1);
for i=1:n
    L(i)=sqrt(sum(pnts(:,i).^2+pnts(:,i+1).^2));
end
q=0;
for i=p+1:n+1
    for j=i-p:i-1
        q=q+L(j);
    end
end
k=n+p+2;
u=zeros(1,k);
for i=p+1:n+1
    s=0;
    for j=i-p:i-1
        s=s+L(j);
    end
    u(i+1)=u(i)+s/q;
end
u(k-p:k)=1;

crv =nrb2nrl( nrbmak(pnts,u) );


%% Demo
% p=2;
% x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
% y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
% crv=nrlspline(p, x, y);
% figure; nrlctrlplot(crv);




