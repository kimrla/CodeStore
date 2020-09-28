function [srfs, crvs]=CreatCoons(ncrvs)

% Create NURL Coons surfaces from NURL or NURBS curves
% 
%   Input:
%      ncrvs - NURL or NURBS curves
%
%  Output:
%      srfs - NURL surfaces
%      crvs -  NURL Coons curves
% 

for i=1:numel(ncrvs)
    if (~isstruct(ncrvs(i)))
        error('NURBS or NURL representation is not structure!');
    end

    if (strcmp(ncrvs(i).form,'B-NURBS'))
        crvs(i)=nrb2nrl(ncrvs(i));
    elseif (strcmp(ncrvs(i).form,'L-NURL'))
        crvs(i)=ncrvs(i);
    else
        error('Not a recognised NURBS or NURL representation');
    end    
end

% Get edges for coons surfaces
coons=edge2coons(crvs);

% Get the coons surfaces
[m,~]=size(coons);
for i=1:m
    tcrvs={crvs(coons(i,1)), crvs(coons(i,2)), crvs(coons(i,3)), crvs(coons(i,4))};
    srfs(i)=nrlcoons(tcrvs{1}, tcrvs{2}, tcrvs{3}, tcrvs{4});
end
srfs=nrlsrfsdirection(srfs);

%% Demo - NURL curves
% load brim sdata   % torquearm
% crvs=sdata.crvs;       
% 
% % Create NURL Coons surfaces from NURL or NURBS curves
% [srfs, crvs]=CreatCoons(crvs);
% 
% % Plot
% figure; hold on;
% for i=1:numel(srfs)
% 	nrlplot(srfs(i),[15,15]);
% end
% chs= nrlcrvplot(crvs, 50);
% shading interp;
% axis equal; 

%% Demo - NURBS curves
% load torquearm sdata   % torquearm
% crvs=sdata.crvs;       
% 
% % Create NURL Coons surfaces from NURL or NURBS curves
% [srfs, crvs]=CreatCoons(crvs);
% 
% % Plot
% figure; hold on;
% for i=1:numel(srfs)
% 	nrlplot(srfs(i),[15,15]);
% end
% chs= nrlcrvplot(crvs, 50);
% shading interp;
% axis equal; 

