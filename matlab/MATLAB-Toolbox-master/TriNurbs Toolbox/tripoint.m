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
% Discription:
%
%     See also:  interline2tri
%

P1=P1(:)'; P2=P2(:)'; P3=P3(:)'; 
w=1-u-v;
dim=length(P1);
if dim==3
    x=P1(1)*u+P2(1)*v+P3(1)*w;
    y=P1(2)*u+P2(2)*v+P3(2)*w;
    z=P1(3)*u+P2(3)*v+P3(3)*w;
    P=[x(:), y(:), z(:)];
elseif dim==2
    x=P1(1)*u+P2(1)*v+P3(1)*w;
    y=P1(2)*u+P2(2)*v+P3(2)*w;
    P=[x(:), y(:)];
end
if nargout>1
    Pu=P1-P3;
    Pv=P2-P3;
    Jac={Pu, Pv};
end





