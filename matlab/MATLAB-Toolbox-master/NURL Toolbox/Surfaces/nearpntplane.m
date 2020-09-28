function [uv, pts, dn, jac]=nearpntplane(srf, pnt)

% Get the nearest point of a plane surface to a given piont
% 
% Calling Sequences:
%
%     [uv, pts, dn, jac]=nearpntplane(srf, pnt)
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

% Get the nearest distance by iteration
[m, n]=size(pnt); m=max([m,n]);
if m>1
    pnt0=pnt; pnt=zeros(1,3);
    pnt(1:m)=pnt0(:);
end
F=zeros(2,1); uv=[0; 0]; dF=zeros(2);
[pts, jac]=nrldeval(srf, uv);
F(1)=sum((pts-pnt').*jac{1});
F(2)=sum((pts-pnt').*jac{2});
dF(1,1)=sum(jac{1}.^2);
dF(1,2)=sum(jac{1}.*jac{2});
dF(2,1)=dF(1,2);
dF(2,2)=sum(jac{2}.^2);
uv=uv-dF\F;
dn=norm(pnt'-pts);


%% ! Demo
% srf = nrl4surf([1.0 0.0 -0.0], [0.0 0.0 0.0], [0.0 1.0 -0.0], [1.0 1.0 0.0]);
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;
% 
% pnt=[0.5, 0.5, 0.5]; hold on;
% plot3(pnt(1), pnt(2), pnt(3), 'ro');
% 
% [uv, pts, dn, jac]=nearpntplane(srf, pnt');
% 
% d=[pnt', pts];
% plot3(d(1,:), d(2,:), d(3,:));
% plot3(pts(1), pts(2), pts(3), 'ro');
% quiver3(pts(1), pts(2), pts(3), jac{1}(1), jac{1}(2), jac{1}(3));
% quiver3(pts(1), pts(2), pts(3), jac{2}(1), jac{2}(2), jac{2}(3));
% dt=cross(jac{1}, jac{2});
% quiver3(pts(1), pts(2), pts(3), dt(1), dt(2), dt(3));






