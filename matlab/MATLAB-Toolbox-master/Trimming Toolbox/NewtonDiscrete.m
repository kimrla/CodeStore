function [pts2,pnts2,d]=NewtonDiscrete(f,pts1,varargin)
% Use Newton-Discrete method, meaning use difference quotient to replace
% the derivative or the element in the Jacabo matrix. This function can use
% Mountain method optionally, and issuitable for both equations and equation.
% The step h=c*||F||, meaning this method is Newton-Steffensen.

% Input:
%   f: Original function handle.
%   pts1, pnts1: Par/phy-coords of the given point.
%   'iterative':  Iimes of iterative, default is 3.
%   'tolerance'=TOL: Tolerance between x(k+1) and x(k), default is eps.
%   'step': Newton-Steffensen step length factor, default is 1e-2
%   'mountain'=MIU: Mountain factor, default is [], meaning NOT use
%       Mountain method.

% Output:
%   pts2: Par-coords of the original function f(x)
%   pnts2: Phy-coords of the original function f(x)
%   d: Distance between the 2 points.

% Default value of the optional variance.
ITE=3;
MIU=[];
TOL=eps;
c=1e-2;% Newton-Steffensen factor, related to the step length.
% c=rand(1,2);

options=struct('iterative',ITE,...
    'mountain',MIU,...
    'tolerance',TOL,...
    'step',c);
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
x0=x0';
pnts1=f(x0); % Output is column vector
for i=1:options.iterative % Endding condition 1.
    dF=discretederiv(f,x0);   
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