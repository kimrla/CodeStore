function [pts2,pnts2,d]=Newton(f,df,pts1,varargin)
% Traditional Newton-Raphson iterative method.This function can use
% Mountain method optionally, and is suitable for both equations and equation..
% Three endding conditions: iterative, tolerance, mountain.

% Input: 
%   f: Original function handle.
%   df: Derivative of original function handle.
%   pts1, pnts1: Par/phy-coords of the given point.
%   'iterative':  Iimes of iterative, default is 5.
%   'mountain'=MIU: Mountain factor, default is [], meaning NOT use
%       Mountain method.

% Output:
%   pts2: Par-coords of the original function f(x)
%   pnts2: Phy-coords of the original function f(x)
%   d: Distance between the 2 points.

% Default value of the optional variance.
ITE=5;
MIU=[];
TOL=eps;

options=struct('iterative',ITE,...
    'mountain',MIU,...
    'tolerance',TOL);
optionnames=fieldnames(options);

nargin=length(varargin);
if round(nargin/2)~=nargin/2
    error('The input is error!');
end

for pair=reshape(varargin,2,[])
    if any(strcmp(pair{1},optionnames))
        options.(pair{1})=pair{2};
    end
end

% Both suitable for equation and equations.
x0=pts1(:);
pnts1=f(x0);
for i=1:options.iterative % Endding condition 1.
    F=f(x0);
    dF=df(x0);  
    % Ending condition 2.    
    if ~isempty(options.mountain)
        x=itemount(f,x0,F,dF,options.mountain);      
    else
        x=x0-dF\F; 
    end     
    % Ending conditionn 3.
    if detercond('tolerance',x0,x,options.tolerance)
        break;   
    else
        x0=x;
    end        
end  

pts2=x;
pnts2=f(x);
d=norm(pnts2-pnts1);

end