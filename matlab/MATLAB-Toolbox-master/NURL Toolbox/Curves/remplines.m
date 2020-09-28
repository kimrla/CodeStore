function [crvs, T]=remplines(crvs)

% Remove lines that are actually points
% 
% Inputs: 
% 
%    crvs - the curves to be used
%
% Output: 
%
%   crvs - the curves without point lines
%   T - indexes of the lines removed
%

% Get a tolerance
n=numel(crvs);
dist=0; 
for i=1:min([10, n])
    disti= nrlmeasure (crvs(i));
    if disti>dist
        dist=disti;
    end
end

% Find the lines that are points
T=ones(1,n);
for i=1:n
    pnts = nrleval(crvs(i), linspace(0, 1, 5)); 
    dm = DistanceMatrix(pnts',pnts'); 
    md=max(max(dm));
    if md<=dist*1e-3
        T(i)=0;
    end
end
L=T==1;
T=find(T==0);

% Remove points
crvs=crvs(L);





