function [Tx,Ty,Txy,Trr,Ttt,Trt,U,V,Ur]=ThickCylinderAnalytical(x,y,p1,p2,a,b,E,nu)

% ThickCylinderAnalytical: Analytical solution for thick-walled cylinder
% 
% Calling Sequences:
%
%     [Tx,Ty,Txy,Trr,Ttt,Trt,U,V,Ur]=ThickCylinderAnalytical(x, y, p1, p2)
%
% INPUTS:
%
%      x, y - Cartesian coordinates of the thick-walled cylinder
%
%      p1, p2 - stress on inner and outter circles
%
%      a, b - inner and outter radius
% 
%      E, nu -  Young's modulus (E) and Poisson's ratio (nu)
%
% OUTPUT:
%
%   Tx, Ty, Txy  -  stress components in Cartesian coordinates
%   Trr, Ttt, Try - stress components in polar coordinates
%   U, V  - displacement components in Cartesian coordinates
%   Ur - displacement components in polar coordinates
% 

if nargin<=4
    a=0.5; b=1;
end
if nargin<=6
    nu=0.3; E=72e9; 
end

r=sqrt(x.^2+y.^2);
t=asin(y./r);

B=a^2*b^2*(p2-p1)/(b^2-a^2);
C=(a^2*p1-b^2*p2)/(b^2-a^2);
Trr=B./r.^2 + C;
Ttt=-B./r.^2 + C;
Trt=0*t;

Ur=((1+nu)/E)*( -(a^2*b^2*(p2-p1)/(b^2-a^2))./r + (1-2*nu)*((a^2*p1-b^2*p2)/(b^2-a^2))*r );

Tx=Trr.*cos(t).^2+Ttt.*sin(t).^2-2*Trt.*cos(t).*sin(t);
Ty=Trr.*sin(t).^2+Ttt.*cos(t).^2+2*Trt.*cos(t).*sin(t);
Txy=(Trr-Ttt).*cos(t).*sin(t)+Trt.*(cos(t).^2-sin(t).^2);

Ut=0*t;
U=Ur.*cos(t)-Ut.*sin(t); 
V=Ur.*sin(t)+Ut.*cos(t);




