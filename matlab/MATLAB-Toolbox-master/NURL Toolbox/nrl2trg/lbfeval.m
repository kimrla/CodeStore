function pp=lbfeval(bsrf, tt)

% Evaluate an lbf surface
% 
% Calling Sequences:
%
%     pp=lbfeval(bsrf, tt)
%
% INPUTS:
%
%      bsrf    :   a blending function triangle.
%
%      tt     - a cell {tu, tv} of parametric evaluation points. 
%                or scattered parametric coordinates (see nrleval)
%
% OUTPUT:
% 
%     pp    :   evaluated points corresponding to the 
%                 parametric points
% 

% Evaluate the surface formed by the edges
[u, v]=meshgrid(tt{1}, tt{2});
S=u.*(1-v); T=v;  
crvs=bsrf.edges;
P1=bsrf.corners(:,1);
P2=bsrf.corners(:,2);
ct=trianglcoef(S(:)', T(:)');
p1=nrleval(crvs(1), S(:)');
p2=nrleval(crvs(2), T(:)');
p3=nrleval(crvs(3), T(:)');
x=ct{1}.*p1(1,:)+ct{2}.*p2(1,:)+ct{3}.*p3(1,:)+ct{4}*P1(1)+ct{5}*P2(1);
y=ct{1}.*p1(2,:)+ct{2}.*p2(2,:)+ct{3}.*p3(2,:)+ct{4}*P1(2)+ct{5}*P2(2);
z=ct{1}.*p1(3,:)+ct{2}.*p2(3,:)+ct{3}.*p3(3,:)+ct{4}*P1(3)+ct{5}*P2(3);
pb=[x; y; z];

% Evaluate the inners error surface
pnts = nrleval(bsrf.efaces, tt);
pnts=permute(pnts, [1 3 2]);
pa=reshape(pnts, 3, []);

% The lbf surface is a summation of the two surfaces above
pp=pb+pa;




