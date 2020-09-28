function [x, pts1, pts2, d, jac1, jac2]=srfcrvintersect(srf, crv, x, it)

% Get intersections of a curve with a surface by Newton-Raphson method
% 
% Calling Sequences:
% 
%     [x, pts1, pts2, d]=srfcrvintersects(srf, crv, x)
% 
% INPUTS:
%
%      srf - a nurls surface
%      crv - a nurls curve
%      x - initial value of parametric intersection point
%      it - the number of iterations
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     nts1 - points of intersections on the surface
%     pts2 - points of intersections on the curve
%     d - distance of the two points
% 

if nargin==3
    it=8;
end

F=zeros(3,1); dF=zeros(3);
for i=1:it
    [pts1, jac1]=nrldeval(srf, x(1:2)'); 
    [pts2, jac2]=nrldeval(crv, x(3));  
    dr=pts1-pts2; 
    F(1)=dot(dr, jac1{1}); 
    F(2)=dot(dr, jac1{2}); 
    F(3)=-dot(dr, jac2); 
    dF(1,1)=dot(jac1{1}, jac1{1}); 
    dF(1,2)=dot(jac1{1}, jac1{2}); 
    dF(1,3)=-dot(jac1{1}, jac2); 
    dF(2,1)=dF(1,2); dF(3,1)=dF(1,3); 
    dF(2,2)=dot(jac1{2}, jac1{2}); 
    dF(2,3)=-dot(jac1{2}, jac2); 
    dF(3,2)=dF(2,3); 
    dF(3,3)=dot(jac2, jac2); 
    x=x-dF\F; 
    pp=x>1; x(pp)=1; 
    pp=x<0; x(pp)=0; 
end
d=norm(pts2-pts1);










