function [uv, pts, dn, jac]=nearpntsrf(srf, pnt, p, q)

% Get the nearest point of a surface to a given piont
% 
% Calling Sequences:
%
%     [uv, pts, dn, jac]=nearpntsrf(srf, pnt)
%     [uv, pts, dn, jac]=nearpntsrf(srf, pnt, p, q)
%
% INPUTS:
%
%      srf - a nurls surface
%      pnt - an arbitrary point
%      p, q - the number of sampling points used get 
%               initial guess of the nearest point
%
% OUTPUT:
% 
%     uv - parametric coordiantes of the nearest point of the surface
% 
%     pts - coordinates of the nearest point of the surface
%
%     dn - the nearest distance from the point tot he surface
%
%     jac - tangent vectors of the point on the surface
% 

if nargin==2
    p=15; q=p;
end

% Get the approximate nearest distance of a point to a surface
[m, n]=size(pnt); m=max([m,n]);
if m>1
    pnt0=pnt; pnt=zeros(1,3);
    pnt(1:m)=pnt0(:);
end
t1=linspace(0, 1, p); 
t2=linspace(0, 1, q); 
pnts=nrleval(srf, {t1, t2});
dm=DistanceMatrix(pnt, pnts(:,:)');
dm=reshape(dm, p, q);
[di, I]=min(dm); 
[~, I2]=min(di);
u=t1(I(I2)); v=t2(I2);

% Get the nearest distance by iteration
[uv, pts, dn, jac]=distpntsrf(srf, pnt, u, v);


%% ! Demo
% R=1; N=6;
% s1=0; s2=pi/2; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;
% 
% pnt=[1, 1, 0.5]; hold on;
% plot3(pnt(1), pnt(2), pnt(3), 'ro');
% 
% [uv, pts, dn, jac]=nearpntsrf(srf, pnt);
% 
% d=[pnt', pts];
% plot3(d(1,:), d(2,:), d(3,:));
% plot3(pts(1), pts(2), pts(3), 'ro');
% quiver3(pts(1), pts(2), pts(3), jac{1}(1), jac{1}(2), jac{1}(3));
% quiver3(pts(1), pts(2), pts(3), jac{2}(1), jac{2}(2), jac{2}(3));
% dt=-cross(jac{1}, jac{2});
% quiver3(pts(1), pts(2), pts(3), dt(1), dt(2), dt(3));







