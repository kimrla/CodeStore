function x=itemount(f,x0,F,dF,varargin)
% The mountatin-down method of iterative.

% Input:
%   f: Original function handle.
%   x0: The (k)th par-coord.
%   pnts1: Phy-coords of the given point.
%   F: Dependent variable f(x(k)).
%   dF: Derivative or the Jacabo matrix of dependent variable df(x(k)).
%   MIU: Mountain factor tolerance, default is 1/2e10.

% Output:
%   x: The (k+1)th par-coord after using Mountain method.

miu=1;
MIU=1/2E10;
if length(varargin)==1
    MIU=varargin{1};
elseif length(varargin)==2
    if strcmp(varargin{1},'mountain')
        MIU=varargin{2};
    else
        error('The input is error!');
    end
end

dr=f(x0);
while ~detercond('mountain',miu,MIU)
    x=x0-miu*dF\F; 
    % Ensure the par-coords are located in the define domain.
    pp=x>1; x(pp)=1; 
    pp=x<0; x(pp)=0;               
    dr2=f(x);
    if norm(dr2)<=norm(dr) % Not larger instead of smaller.
        break;
    else
        miu=miu/2;                 
    end              
end

end