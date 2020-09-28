function [pts2, pnts2, dist]=nrbsrfreversebi(srf, x, pnts1, it, n)

% tnrbintersect: Get the intersection points of two tri-nurbs surfaces by method of bisection
% 
% Calling Sequences:
% 
%       [pts2, pnts2, dist]=nrbsrfreversebi(srf, x, pnts1)
% 
%       [pts2, pnts2, dist]=nrbsrfreversebi(srf, x, pnts1, it)
% 
% INPUTS:
% 
%       srf - A nurbs surface.
%
%       x=[u, v]   -  Initial value of the parametric coordinates.
% 
%       pnts1 - A given point on Cartesian coordinate.
%
%       it -  The number of iterations. The default value is 3.
%
%       n - The number of seeds.
% 
% OUTPUT: 
% 
%       pts2 - Parametric points of the surface.
% 
%       pnts2 - The nearest points of the surface to given points.
% 
%       dist  -  Distances between points on the surface and the given points.
%

if nargin==3
    it=10; n=8;
end
if nargin==4
    n=8;
end

ds=0.5; 
x=x(:); pnts1=pnts1(:)';
for i=1:it
    s=linspace(max([0, x(1)-ds]), min([1, x(1)+ds]), n);
    t=linspace(max([0, x(2)-ds]), min([1, x(2)+ds]), n);
    pt2=nrbeval(srf, {s, t});
    pt2=pt2(:,:)'; 
    [~, pti]=nearpnts(pnts1, pt2);
    pnts2=pt2(pti, :); 
    i2=rem(pti, n); 
    j2=fix(pti/n)+1; 
    if i2==0
        j2=j2-1; i2=n;
    end
    x=[s(i2), t(j2)];
    ds=ds/2;
end
pts2=x;
dist=norm(pnts2-pnts1);


%% demo
% % Create a nurbs surface and a point
% circ=nrbcirc(4, [5,5,4], 0, pi);
% srf1=nrbrevolve(circ, [5,5,4], [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% pts=nrbeval(srf1, {[0.2, 0.4, 0.8]; [0.2, 0.4, 0.8]});
% pnts1=pts(:,:)';
% srf=srf1;
% 
% % The nearest points of a surface to given points
% s=linspace(0,1,10); t=linspace(0,1,11);
% [t,s]=meshgrid(t,s);
% pts=[s(:)'; t(:)'];
% points=nrbeval(srf, pts)';
% [~, p2]=nearpnts(pnts1, points);
% np=length(pnts1);
% pnts2=zeros(np,3);
% pts2=zeros(np,2);
% dist=zeros(np,1);
% dsrf=nrbderiv(srf);
% for k=1:np
%     [pts2(k,:), pnts2(k,:), dist(k)]=nrbsrfreversebi(srf, pts(:, p2(k)), pnts1(k,:)');
% end
% 
% % Plot the results
% figure; hold on;
% nrbplot(srf, [20, 20]);
% plot3(pnts1(:,1), pnts1(:,2), pnts1(:,3), 'ro');
% plot3(pnts2(:,1), pnts2(:,2), pnts2(:,3), 'r*');
% axis equal; view(3);
% title('Geometric grid');





