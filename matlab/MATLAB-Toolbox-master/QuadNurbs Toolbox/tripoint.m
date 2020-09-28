 function [P, Jac]=tripoint(P1, P2, P3, u, v)

% tripoint: Get points on a triangle.
% 
% Calling Sequences:
% 
%     P=tripoint(P1, P2, P3, u, v)
% 
%     [P, Jac]=tripoint(P1, P2, P3, u, v)
% 
% INPUTS:
%
%      P1, P2, P3 - Coordinates of the three vertexes of the triangle.
%
%      u, v - Parametric coordinates of the triangle defined on [0, 1].
%
% OUTPUT:
% 
%     P - Points corresponds to the parametric coodinates.
%
%     Jac = [Pu, Pv]  -  Directions on the two parametric coodinates.
%

P1=P1(:)'; P2=P2(:)'; P3=P3(:)'; 
w=1-u-v;
dim=length(P1);
% Use area coordinates of the triangle to obtain as the local coordinates
if dim==3
    x=P1(1)*w+P2(1)*u+P3(1)*v;
    y=P1(2)*w+P2(2)*u+P3(2)*v;
    z=P1(3)*w+P2(3)*u+P3(3)*v;
    P=[x(:), y(:), z(:)];
elseif dim==2
    x=P1(1)*w+P2(1)*u+P3(1)*v;
    y=P1(2)*w+P2(2)*u+P3(2)*v;
    P=[x(:), y(:)];
end
if nargout>1
    Pu=P2-P1;
    Pv=P3-P1;
    Jac={Pu, Pv};
end





