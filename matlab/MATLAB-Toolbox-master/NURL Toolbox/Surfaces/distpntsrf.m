function [uv, pts, dn, jac]=distpntsrf(srf, pnt, u, v, it)

% Get the nearest point of a surface to a given piont
% 
% Calling Sequences:
%
%     [uv, pts, dn, jac]=distpntsrf(srf, pnt, u, v)
%     [uv, pts, dn, jac]=distpntsrf(srf, pnt, u, v, it)
%
% INPUTS:
%
%      srf - a nurls surface
%      pnt - an arbitrary point
%      u, v - approximate parametric coordiantes of 
%                the nearest point of the surface
%      it - the number of iterations
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
%  See also:  nearpntsrf
%

if nargin==4
    it=5;
end

% Get the nearest distance by iteration
F=zeros(2,1); uv=[u; v]; dF=zeros(2);
for i=1:it
    [pts, jac, hess]=nrldeval(srf, uv);
    F(1)=sum((pts-pnt').*jac{1});
    F(2)=sum((pts-pnt').*jac{2});
    dF(1,1)=sum(jac{1}.^2)+sum((pts-pnt').*hess{1,1});
    dF(1,2)=sum(jac{1}.*jac{2})+sum((pts-pnt').*hess{1,2});
    dF(2,1)=dF(1,2);
    dF(2,2)=sum(jac{2}.^2)+sum((pts-pnt').*hess{2,2});
    uv=uv-dF\F;
end
dn=norm(pnt'-pts);






