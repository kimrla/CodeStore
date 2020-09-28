function [crvs1, crvs2, uu, pp, d]=intsctpnts(crv1, crv2, n)

% Solve points of intersection of two curves
% 
%  Inputs: 
%
%    crv1, crv2 - the two curves to be used
% 
%    n - the number of nodes (the default value is 100)
%
%  Output: 
%     crvs1 - curves splitted from crv1
%     crvs2 - curves splitted from crv2
%     uu - the intersection parametric points
%     pp - the coordinates of the intersection points of the two curves
%     d - the minimum distances of the two curves
% 
%  Examples:
%     p=2; pnt=[6;5.4;0]; 
%     crv1=nrlcirc(3, pnt', pi, 2*pi);
%     x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
%     y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
%     crv2=nrlspline(p,x,y); 
%     [crvs1, crvs2, uu, pp, d, flag]=intsctpnts(crv1, crv2);
%     [crvs1, crvs2, uu, pp, d, flag]=intsctpnts(crv1, crv2, 100);
% 

if nargin==2
    n=100;
end

% Get the length of the curves
dist= max([nrlmeasure(crv1), nrlmeasure(crv2)]);

% Get approximated points of intersection of two curves
[uu, du] = aintsctpnts(crv1, crv2, n);

% Solve points of intersection using Newton-Raphson's method
m=length(du); d=zeros(m,1); pp=zeros(3,2,m); 
for i=1:m
    [uu(i,:), d(i), pp(:,:,i)]=intsctpnt(crv1, crv2, uu(i,:)');
end
id=d/dist<1e-6;
d=d(id); uu=uu(id,:); pp=pp(:,:,id);
if isempty(d)
    crvs1=crv1;
    crvs2=crv2;
    return;
end

% Split the curves
crvs1=nrlsplits(crv1, uu(:,1));
crvs2=nrlsplits(crv2, uu(:,2));
crvs1=remplines(crvs1);
crvs2=remplines(crvs2);


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
% [crvs1, crvs2, uu, pp, d]=intsctpnts(crv1, crv2, n);
% 
% figure; hold on;
% for i=1:numel(crvs1)
%     nrlplot(crvs1(i));
% end
% for i=1:numel(crvs2)
%     nrlplot(crvs2(i));
% end
% 
% for i=1:length(d)
%     plot(pp(1,1,i), pp(2,1,i), 'ro');
%     plot(pp(1,2,i), pp(2,2,i), 'bs');
% end


