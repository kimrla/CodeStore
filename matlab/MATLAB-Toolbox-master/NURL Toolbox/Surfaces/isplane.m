function isp=isplane(srf)

% Check whether a surface is a plane
% 
% Calling Sequences:
%
%     isp=isplane(srf)
%
% INPUTS:
%
%      srf - a nurls surface
%
% OUTPUT:
% 
%     isp - true (1) or false (0)
% 

p=5; q=6;
t1=linspace(0, 1, p); 
t2=linspace(0, 1, q); 
[~, jac]=nrldeval(srf, {t1, t2});
dn=cross(jac{1}(:,:), jac{2}(:,:));
dn=vecnorm(dn);
dm = DistanceMatrix(dn', dn');
di=max(max(dm));
if di>1e-10
    isp=0;
else
    isp=1;
end


%% ! Demo
% srf = nrl4surf([1.0 0.0 -0.0], [0.0 0.0 0.0], [0.0 1.0 -0.0], [1.0 1.0 0.0]);
% isp=isplane(srf);


