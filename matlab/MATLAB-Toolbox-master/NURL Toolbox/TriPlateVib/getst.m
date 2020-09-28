function [S, T, u, v]=getst(tt, di)

% Get the matrices of u and v coordinates for the triangle
% 
% Calling Sequences:
%
%     [S, T]=getstc(tt, di)
%
% INPUTS:
% 
%      tt     - parametric evaluation points. The nurls should be  
%            a surface and tt is a cell {tu, tv} of the parametric coordinates
%
%     di -  direction of integration 1 (u) or 2 (v)
%
% OUTPUT:
% 
%     S    :   a matrix of u coordinates for integration on a triangle
%     T    :   a matrix of v coordinates for integration on a triangle
%     u, v :  matrices of u,v coordinates for on a rectangular domain
% 

if nargin==1
    di=1;
end

if iscell(tt)
    [v, u]=meshgrid(tt{2}, tt{1});
else
    [m,~]=size(tt);
    if m~=2
        tt=tt';
    end
    u=tt(1,:);
    v=tt(2,:);
end

switch di
    case 1        
        S=u.*(1-v);
        T=v;        
    case 2    
        S=u;
        T=v.*(1-u);
end





