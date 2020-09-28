function [x, pnts1, pnts2, dt]=optintersects(x, pnts1, pnts2, dt, tol)

% Remove duplicated intersections and sort the intersection points
% 
% Calling Sequences:
% 
%     [x, pnts1, pnts2, dt]=optintersects(x, pnts1, pnts2, dt)
% 
% INPUTS:
%
%     x - parametric points of intersections
%     pnts1 - points of intersections on surface 1
%     pnts2 - points of intersections on surface 2
%     dt - distance of the two nearest points
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     pnts1 - points of intersections on surface 1
%     pnts2 - points of intersections on surface 2
%     dt - distance of the two nearest points
% 
%  See also:
%       srfsinterscts
%

if nargin==4
    tol=1e-3;
end

% Remove duplicated intersections
[x, uu]=RemDuplicate(x', tol);
x=x';
pnts1=pnts1(:,uu);
pnts2=pnts2(:,uu);
dt=dt(uu);

% Sort the intersection points
n=length(dt);
[~, I]=min(min(x,[],2));
if I==1 || I==2
    [~, p]=sort(x(I,:)); ps=zeros(size(p)); 
    dm=DistanceMatrix(x(1:2, p)', x(1:2, p)'); 
else
    [~, p]=sort(x(I,:)); ps=zeros(size(p)); 
    dm=DistanceMatrix(x(3:4, p)', x(3:4, p)'); 
end
ps(1)=p(1); 
id=true(n,1);
ip0=1;
for i=2:n
    id(ip0)=false;
    dd=dm(ip0, id);
    [~, ip] = min(dd);
    pp=find(id); 
    ip=pp(ip);
    ps(i)=p(ip); 
    ip0=ip;
end
x=x(:,ps);
pnts1=pnts1(:,ps);
pnts2=pnts2(:,ps);
dt=dt(ps);
ps=find(x(I,:)>=1, 1, 'first')+1;
pnts1(:,ps:end)=[];
pnts2(:,ps:end)=[];
x(:,ps:end)=[];
dt(ps:end)=[];






