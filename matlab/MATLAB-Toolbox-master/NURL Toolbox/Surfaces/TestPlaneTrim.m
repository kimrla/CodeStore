clear; clc;

% Get a plane and insert knots
a=5;
srf = nrl4surf([a 0.0 0], [0.0 0.0 0], [0.0 2.0 0], [a 2.0 0]);
iuknt=20*ones(1, length(srf.intervals{1})-1);
ivknt=10*ones(1, length(srf.intervals{2})-1);
srf = nrlkntins(srf, {iuknt, ivknt});
figure; nrlctrlplot(srf);
axis equal; view(2); hold on;

% Plot curves on the plane
crv1=nrlcirc(1, [a-1, 1], -pi/2, pi/2);
crv2=nrlcirc(1, [1, 1], pi/2, 3*pi/2);
crv3=nrlcirc(0.5, [a-1, 1]);
crv4=nrlcirc(0.5, [1, 1]);
crv5=nrlline([1,0], [a,0]);
crv6=nrlline([1,2], [a-1,2]);
crvs=[crv1, crv2, crv3, crv4, crv5, crv6];
nrlcrvplot(crvs, 100, 'k'); 

% Get the parametric points on the plane surface
pps=cell(5,1); m=35;
t=linspace(0, 1, m);
uv=zeros(2,m);
for i=1:numel(crvs)
    pts=nrleval(crvs(i), t);
    for j=1:m
        uv(:,j)=nearpntplane(srf, pts(:,j));
    end
    pps{i}=uv;
end

% Plot parametric points
figure; hold on;
[u, v]=meshgrid(srf.knots{1}, srf.knots{2});
plot(u, v); plot(u', v');
for i=1:numel(crvs)
    plot(pps{i}(1,:), pps{i}(2,:));
end
axis equal;



