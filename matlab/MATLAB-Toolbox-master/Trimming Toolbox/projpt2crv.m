function [pts2,pnts2,d]=projpt2crv(crv,pts1,pnts1,varargin)
% Projection from one given point to the NURBS curve, or inverse the parametric
% coords of the given physical coorde, using Newton-Mountain method,
% comparing with projpt2srf.m. Just consider the normal condition: The
% derivative can be obtained directly.

% Three endding conditions: iterative, tolerance, mountain.

% Input: 
%   crv: NURBS/Bezier structure curve.
%   pts1, pnts1: Par/phy-coords of the given point.
%   'iterative':  Iimes of iterative, default is 3.
%   'tolerance'=TOL: Tolerance between x(k+1) and x(k), default is eps.
%   'mountain'=MIU: Mountain factor, default is 1/2e10

% Output:
%   pts2: Par-coords of the prjoection point in the surface.
%   pnts2: Phy-coords of the prjoection point in the surface.
%   d: Distance between the 2 points.

% Default value of the optional variance.
ITE=3;
MIU=1/2e10;
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

dcrv=nrbderiv(crv);   
% dr=nrbeval(crv,knt)-pnts1, size(dr)=[3,1]
% Handle of original function, size(f(x))=[1,1], size(df(x))=[1,1]
f=@(knt) (jacobi(crv,dcrv,knt)*(nrbeval(crv,knt)-pnts1));
df=@(knt) (jacobi(crv,dcrv,knt)*jacobi(crv,dcrv,knt)');
[pts2,pnts2,d]=Newton(f,df,pts1,'iterative',ITE,...
                                    'mountain',MIU,...
                                    'tolerance',TOL); 
end

function dF=jacobi(crv,dcrv,knt)

[~,dF]=nrbdeval(crv,dcrv,knt); 
if iscell(dF)
    dF=[dF{:}];% size dF=[3,2]
end
dF=dF'; % size(dF)=[2,3] or [1,3]

end




