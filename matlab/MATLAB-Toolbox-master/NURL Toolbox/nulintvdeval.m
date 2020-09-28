function dpnts=nulintvdeval(t, knots, coefs, order, der)

% Evaluate nul curves or their derivatives at an interval 
% 
%   input : 
%               t - the evalued point 
%               knots - knots vector
%               coefs - coefs matrix
%               order - the order of the Langrange basis
%               der - the order of derivatives
% 
%  output :   
%    dpnts : Evaluated points on the NUL curve as Cartesian coordinates (x,y,z)
%

% Check order
if order>7
    order=7;
end

if isempty(t)
    dpnts=[];
    return;
end

% Get span index
[index, iu]=spanweak(knots, t, order);

% Get the points
m=length(t); [mc, ~] = size(coefs);
dpnts=zeros(mc, m); 
[n, ~]=size(iu); 
for i=1:n
    pp=index(i)+1:index(i+1);
    dpnts(:, pp)=nulspandeval(t(pp), coefs(:,iu(i,:)), knots(iu(i,:)), order, der); 
end


function pnts=nulspandeval(x, coefs, knots, order, der)

% Evaluate nul curve or their derivatives at a span
% 
%   input :  x - the evalued point
%               knots - knot vector of each basis
%               order - the order of the Langrange basis
%               der - the order of derivatives
% 
%  output :   
%    pnts : Evaluated points on the NUL curve as Cartesian coordinates (x,y,z)
%

% Check order
if order>7
    order=7;
end

% Get the length of knots (n) and the number of 
%   points (m) to be evaluated
n=length(knots); m=length(x);
[mc, ~] = size(coefs);
pnts=zeros(mc, m);

% Check whether the order agrees with the length 
%   of the knots
if order ~= n-1
    error('The knots is not correct!');
end

% Get the values of the basis at x
dnL=nurlbasis(x, knots, order, der);

% The values of the points                
for i=1:mc
    v=0;
    for j=1:n
        v=v+coefs(i, j)*dnL(j,:);
    end
    pnts(i,:)=v;
end



