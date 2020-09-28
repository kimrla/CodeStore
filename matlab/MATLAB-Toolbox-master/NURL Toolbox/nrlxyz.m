function [x, y, z] = nrlxyz(nrl, mn)

% NRLXYZ: Get the x, y, z coordinates of a nurl surface.
% 
% INPUT:
%
%    srf :  a nurl surface.
%    mn :  the number of grid points in x- (and y)-direction
%
% OUTPUT:
%
%    x :  a m*n matrix of the x coordinates of the nurl surface. 
%    y :  a m*n matrix of the y coordinates of the nurl surface. 
%    z :  a m*n matrix of the z coordinates of the nurl surface. 
% 
% Description:
%     
%    Please read nrlsrfplot fucntion for the usage of this function.
% 

if length(mn)==2
    m=mn(1); n=mn(2);
    
    s=linspace(0,1,m); t=linspace(0,1,n);
    s=[s, nrl.intervals{1}]; s=RemDuplicate(s'); s=sort(s)';
    t=[t, nrl.intervals{2}]; t=RemDuplicate(t'); t=sort(t)';
    p = nrleval (nrl, {s, t});
    x=zeros(n,m); y=x; z=x;
    for j=1:n
        for i=1:m
            x(j,i)=p(1,i,j);
            y(j,i)=p(2,i,j);
            z(j,i)=p(3,i,j);
        end
    end
elseif length(mn)==1
    m = mn; 
    u = linspace(0, 1, m); 
    u=[u, nrl.intervals]; u=RemDuplicate(u'); u=sort(u)';
    p = nrleval(nrl, u); 
    x = p(1,:); 
    y = p(2,:); 
    z = p(3,:); 
else
    error('The number of grid points is not right.');
end