function x=planelineints(pts1, jac1, pts2, jac2, x)

% Get intersections of a straight line with a plane
% 
% Calling Sequences:
% 
%     [x, pts1, pts2, d]=planelineints(srf, crv, x)
% 
% INPUTS:
%
%      srf - a nurls surface
%      crv - a nurls curve
%      x - initial value of parametric intersection point
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     nts1 - points of intersections on the surface
%     pts2 - points of intersections on the curve
%     d - distance of the two points
% 

F=zeros(3,1); dF=zeros(3);

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










