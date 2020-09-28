function crv=nulpts2crv(pnts, order)

% Get a nul curve through arbitrary points
% 
% Input:
%     pnts - points compose of x, y, z coordinates
%     order - the order of the curve (default value is 2)
% 
% Output:
%     crv - a nul curve

if nargin==1
    order=2; 
end

[~, m]=size(pnts);

% Get the distance between points
d=zeros(m-1,1);
for i=1:m-1
   d(i)=norm(pnts(:,i+1)-pnts(:,i));
end

% Get knot vector
u=zeros(1,m);
t=0; sd=sum(d);
for i=2:m
    t=t+d(i-1);
    u(i)=t/sd;
end

% Make a nul curve
if m==2
    order = 1;
end
op=ones(1,m);
coefs=[pnts; op];
crv=nrlmake(coefs, u, [0, 1], order);


%% Demo
% % Points on a plane curve
% pts=[1, 0.98, 0.92, 0.825, 0.697, 0.54
% 0, 0.199, 0.389, 0.565, 0.717, 0.841
% 0 0 0 0 0 0];
% 
% % Get the curve
%  crv=nulpts2crv(pts);
%  
%  % Get the points
%  n=20;
% t=linspace(0, 1, n);
% pnts=nulintvdeval(t, crv.knots, crv.coefs, crv.order, 0);
% 
% % Plot the curve
% figure; hold on;
% plot(pnts(1,:), pnts(2,:));
% plot(pts(1,:), pts(2,:), 'ro');



