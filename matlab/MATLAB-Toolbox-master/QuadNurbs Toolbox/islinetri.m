function [t, tol]=islinetri(T)

% islinetri: Check whether a triangle is a line
% 
% Calling Sequences:
% 
%     t=islinetri(T)
% 
% INPUTS:
%
%      T=[P1x, P1y, P1z
%           P2x, P2y, P2z
%           P3x, P3y, P3z]
%           Three vertexes of a triangle.
% 
% OUTPUT:
% 
%     t - A logical variable. t=1 if the triangle is a line.
%          Otherwise t=0.
%
%     tol - Maximum edge length of the triangle.
% 

t=false;
tol=max([norm(T(1,:)-T(2,:)), norm(T(1,:)-T(3,:)),norm(T(3,:)-T(2,:))]);
if norm(T(1,:)-T(2,:))<tol*1e-10
    t=true;
end
if norm(T(1,:)-T(3,:))<tol*1e-10
    t=true;
end
if norm(T(3,:)-T(2,:))<tol*1e-10
    t=true;
end






