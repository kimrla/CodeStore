function [pnt, cvt]=curvature(crv, t)

% Get the curvature of a curve
%
% Calling Sequences:
%
%     cvt=curvature(crv, t)
%
% INPUTS:
%
%      crv - a nurl curve
% 
%      t     - parametric evaluation points.
%
% OUTPUT:
%
%   cvt  - curvature of the curve ar parametric points
%

[pnt, jac, hess] = nrldeval (crv, t);
n=length(t);
cvt=zeros(1,n);
for i=1:n
    cvt(i) = sqrt(( sum(jac(:,i).^2)*sum(hess(:,i).^2) - sum(jac(:,i).*hess(:,i))^2 )/( sum(jac(:,i).^2)^3) );
end


%% ! Demo
% % Get elliptic arcs and elevate degrees
% a=2; b=1; N=6;
% sang=0; eang=pi;
% center=[0, 0];
% crv = nrlellip(a, b, center, sang, eang);
% figure; hold on;
% nrlplot(crv, 20);
% axis equal;
% 
% % Get curvature
% t=linspace(0, 1, 11);
% [pnt, cvt]=curvature(crv, t);
%



