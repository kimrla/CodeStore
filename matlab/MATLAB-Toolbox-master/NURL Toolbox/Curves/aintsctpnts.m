function [uu, d, pp] = aintsctpnts(crv1, crv2, n)

% Get approximated points of intersection of two curves
%              using midpoint finite difference method
% 
%  Inputs: 
%
%     crv1, crv2 - the two curves to be used
% 
%     n - the number of nodes (the default value is
%               computed from curvature)
%
% Output: 
%
%     uu - the approximated intersection parametric points
%
%     pp - the approximated coordinates of the intersection 
%                            points of the two curves
% 
%     d - the minimum distances of the two curves
%
% Examples:  
%  
%     p=2; pnt=[6;5.4;0]; 
%     crv1=nrlcirc(3, pnt', pi, 2*pi);
%     x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
%     y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
%     crv2=nrlspline(p,x,y); 
%     [uu, d, pp] = aintsctpnts(crv1, crv2);
%     [uu, d, pp] = aintsctpnts(crv1, crv2, 100);
% 

if nargin==2
    n=100;
end

% Get the points of the curves
ut=linspace(0,1,n);
pnts1=nrleval(crv1, ut);   
pnts2=nrleval(crv2, ut);   

% Get the index of approximated nearest points according to distance
dm = DistanceMatrix(pnts1',pnts2'); 
[midm, II]=min(dm);
dmidm=mfdm(ut, midm);
t=1; IJ=[]; 
for i=1:n-1
    if dmidm(i)<0 && dmidm(i+1)>0
        IJ(t,1)=II(i+1); 
        IJ(t,2)=i+1;
        t=t+1;
    end
end
if isempty(IJ)
    [~, i]=min(midm);
    IJ(t,1)=II(i); 
    IJ(t,2)=i;
end
uu=[ut(IJ(:,1))', ut(IJ(:,2))'];    

% Get the coordinates of approximated nearest points
[m,~]=size(IJ); d=zeros(m,1);
pp=zeros(3,2,m); 
for i=1:m
    d(i)=sqrt(sum((pnts1(:,IJ(i,1))-pnts2(:,IJ(i,2))).^2));
    pp(:,:,i)=[pnts1(:,IJ(i,1)), pnts2(:,IJ(i,2))];
end

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
% n=100;
% [uu, d, pp] = aintsctpnts(crv1, crv2, n);
% 
% for i=1:length(d)
%     plot(pp(1,1,i), pp(2,1,i), 'ro');
%     plot(pp(1,2,i), pp(2,2,i), 'bs');
% end



