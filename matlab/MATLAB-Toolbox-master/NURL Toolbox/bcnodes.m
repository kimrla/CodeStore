function Eg=bcnodes(BC, M, N)
% 
% bcnodes: Get a logical vector of boundary nodes (2D)
% 
% Calling Sequence:
% 
%   Eg=bcnodes(BC, M, N)
% 
% INPUT:
% 
%   BC	: boundary conditions: 1 - clamped, 0 - free
% 
%   center	: Center of the circle, default (0,0,0)
% 
%   M, N	: The number of nodes in x and y diretions
% 
% OUTPUT:
%
%   Eg		:  A logical vector of boundary nodes
%

ex=ones(1, M); ey=ones(1, N); 
gx=0*ex; gx(1)=BC(1); gx(end)=BC(3); 
gy=0*ey; gy(1)=BC(2); gy(end)=BC(4); 
Eg1=kron(ey, gx);
Eg2=kron(gy, ex);
Eg=logical(Eg1) | logical(Eg2);

