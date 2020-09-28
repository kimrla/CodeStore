clear; clc;

load brim sdata   % torquearm
crvs=sdata.crvs;       

% Create NURL Coons surfaces from NURL or NURBS curves
[srfs, crvs]=CreatCoons(crvs);

% Plot
figure; hold on;
for i=1:numel(srfs)
	nrlplot(srfs(i),[15,15]);
end
chs= nrlcrvplot(crvs, 50);
shading interp;
axis equal; 




