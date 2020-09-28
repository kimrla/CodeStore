function [npts,npnts,d]=interbidi(srf1,srf2,pts1,pnts1,pts2,pnts2,varargin)
% Use the binary division method to refine the intersection points.

% Input:
%   srf1, srf2: NURBS/Bezier structure surface.
%   pts1,pnts1: Par/phy-coords of the 1st intersection point.
%   pts2,pnts2: Par/phy-coords of the 2nd intersection point.
%   'mode': 3-Use 3-par iterative method to get the precise point; 
%       4(default)-Use 4-par method.
%   ite: Iterative times, the default equals to 3.
%   method: The default method is Newton-Raphson iterative method And if 
%       input 'mountain', then use Mountain method compulsorily.If the 
%       derivative cannot be calculated (eg. irregular surface), then
%       automatically use Secant method. 
%   st2: If mode='3', [Signal, s or t] of parameter s or t, to determine 
%       which parameter is glued. The first elemtn 1 for gluing the 1st  
%       par(s), 2 for gluing the 2nd par(t). The 2nd element is the glued parameter.

% Output:
%   npts,npnts: New intersection point's PRECISE par/phy-coords by using
%       binary division method to get it.
%   d: Distance between the intersection point npnts's 2 projection points
%       on srf1 and srf2.

ITE=3;
MIU=[];
TOL=eps;
ST2=[]; % Specifically for 3-par iterative method.
MODE=4;

options=struct('iterative',ITE,...
    'mountain',MIU,...
    'tolerance',TOL,...
    'st2',ST2,...
    'mode',MODE);
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

pts=(pts1+pts2)/2;
pnts=(pnts1+pnts2)/2;
if options.mode==3 % 3-par method, ite3par.m
    [npts,npnts,d]=ite3par(srf1,srf2,options.st2,pts,pnts,...
        'iterative',options.iterative,...
        'mountain',options.mountain,...
        'tolerance',options.tolerance);
else % mode=4, 4-par method, ite4par.m
    [npts,npnts,d]=ite4par(srf1,srf2,pts,pnts,...
        'iterative',options.iterative,...
        'mountain',options.mountain,...
        'tolerance',options.tolerance);
end

end