function [pnts1, pnts2, x, d, inter]=interline2intri(T, L)
%
% interline2intri: Get the intersection of a line to the interior of a triangle
% 
% Calling Sequences:
% 
%     [pnts1, pnts2, x, d, inter]=interline2intri(T, L)
% 
% INPUTS:
%
%      T=[P1x, P1y, P1z
%           P2x, P2y, P2z
%           P3x, P3y, P3z]
%           Three vertexes of a triangle.
%
%      L=[P1x, P1y, P1z
%           P2x, P2y, P2z]
%           Two vertexes of a line:
% 
% OUTPUT:
% 
%     pnts1 - Points of intersections on the triangle.
% 
%     pnts2 - Points of intersections on the line.
% 
%     x =[u; v; s], where [u; v] are parametric points of 
%          intersection for the triangle, and s is that of the line.
% 
%     d - Distance of the line to the triangle.
%
%     inter - A value that indicates whether the line is parallel
%              (inter==0) or intersected with the triangle (inter>0). 
%

[t, tol]=islinetri(T);
if t
    pnts1=zeros(0,3);  
    pnts2=zeros(0,3); 
    x=zeros(3,0); 
    d=[]; inter=[];
    return;
end
[pnts1, pnts2, x, d, inter]=interline2tri(T, L);
p=x(3,:)>-eps & x(3,:)<1+eps; 
p=p & d<tol*1e-10;
p=p & x(1,:)>-eps & x(1,:)<1+eps; 
p=p & x(2,:)>-eps & x(2,:)<1+eps; 
p=p & x(1,:)+x(2,:)>-eps & x(1,:)+x(2,:)<1+eps; 
x=x(:,p); 
pnts1=pnts1(p,:); 
pnts2=pnts2(p,:); 






