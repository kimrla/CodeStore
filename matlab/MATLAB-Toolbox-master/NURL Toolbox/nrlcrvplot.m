function fh = nrlcrvplot(crvs, m, clr)

% NRBSRFPLOT: Plot nurl curves and return their figure handles.
% 
% INPUT:
%
%    crvs :  nurl curves.
%    m :  the number of grid points
%
% OUTPUT:
%
%    fh :  handle of the curve
% 


% Draw the surface
if nargin ==2
    clr = 'k';
end

% Draw the curves
n=length(crvs);
fh=zeros(1,n);
for i=1:n
    [x, y, z] = nrlxyz(crvs(i), m);
    fh(i)=line(x, y, z, 'color', clr);
end


    
    