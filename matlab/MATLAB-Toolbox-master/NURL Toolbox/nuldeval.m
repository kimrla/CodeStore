function dpnts=nuldeval(order, cpnts, knts, t, intervals, der)

% Evaluate nul curve, surface, volume or their derivatives
% 
%   input :  
%               order - the order of the Langrange basis
%               cpnts - control points or weighting coefficients, matrix of size (dim,nc). 
%               knts - knots vetor
%               t - the evalued parametric points
%               intervals - intervals of the nuls
%               der - the order of derivatives
% 
%  output :   
%    dpnts : Evaluated points on the NUL curve as Cartesian coordinates (x,y,z)
%                and the corresponding derivatives with respect to
%                parametric points
% 

if nargin==5
    der=0;
end

[mc, ~] = size(cpnts); 
n=length(t); 
dpnts=zeros(mc, n); 
m=length(intervals); 
for i=1:m-1
    a=intervals(i); b=intervals(i+1); 
    p=knts>=a & knts<=b; 
    q=t>=a & t<=b; 
    
    dpntsi=nulintvdeval(t(q), knts(p), cpnts(:,p), order, der); 
    if ~isempty(dpntsi)
        dpnts(:, q)=dpntsi(:, :); 
    end
end

% The case when t<a
a=intervals(1); b=intervals(2); 
p=knts>=a & knts<=b; 
q=t<a; 
dpntsi=nulintvdeval(t(q), knts(p), cpnts(:,p), order, der); 
if ~isempty(dpntsi)
    dpnts(:, q)=dpntsi(:, :); 
end

% The case when t>b
a=intervals(end-1); b=intervals(end); 
p=knts>=a & knts<=b; 
q=t>b; 
dpntsi=nulintvdeval(t(q), knts(p), cpnts(:,p), order, der); 
if ~isempty(dpntsi)
    dpnts(:, q)=dpntsi(:, :); 
end




