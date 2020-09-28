crvs = nrbtestcrv;
crv=crvs;

figure; nrbctrlplot(crvs);

% Transform nurbs curve into nurls curve
crvl=nrb2nrl(crv);

% Evalueate the curves and their derivatives
figure; hold on;
nrlplot(crvl, 1000, 'ctrl'); 
axis equal;

% Load surfaces
load torque srfs
figure; hold on;
for i=1:numel(srfs)
    srfs(i).coefs(3,:,:)=0;
    srfl(i)=nrb2nrl(srfs(i));
    nrlplot(srfl(i), [10, 10],'ctrl');
end
axis equal;

for i=1:numel(srfl)
    nrledgeplot(srfl(i), 10)
end

% Extrude the surfaces
figure; hold on;
for i=1:numel(srfs)
    vols(i)=nrlextrude(srfl(i), [0 0 0.1]);
    nrlplot(vols(i), [10, 10, 10], 'ctrl');
    nrledgeplot(vols(i), 10)
end
view(3); axis equal;