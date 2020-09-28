function [pnts2, pnts1, x]=tnrbinterline2tri(tnrb1, tnrb2, ed1, tri2)


% tnrbinterline2tri: Get the intersection of an edge of tnrb1 to a triangle of tnrb2
% 
% Calling Sequences:
% 
%     [pnts1, pnts2, x]=tnrbinterline2tri(tnrb1, tnrb2, ed1, tri2)
% 
% INPUTS:
% 
%       tnrb1, tnrb2 - Triangular representation of nurbs surface.
%
%      ed1  -  An edge of tnrb1.
%
%      tri2  -  A triangle of of tnrb2.
% 
% OUTPUT:
% 
%     pnts1 - Points of intersections on the triangle.
% 
%     pnts2 - Points of intersections on the line.
% 
%     x =[u, v, s], where [u, v] are parametric points of 
%          intersection for the triangle, and s is that of the line.
%

L=tnrb1.points(ed1, :); 
T=tnrb2.points(tri2, :);
[pnts2, pnts1, x]=interline2tri(T, L);
p=x(3,:)>-eps & x(3,:)<1+eps; 
p=p & x(1,:)>-eps & x(1,:)<1+eps; 
p=p & x(2,:)>-eps & x(2,:)<1+eps; 
p=p & x(1,:)+x(2,:)>-eps & x(1,:)+x(2,:)<1+eps; 
x=x(:,p); 
pnts1=pnts1(p,:); 
pnts2=pnts2(p,:); 






