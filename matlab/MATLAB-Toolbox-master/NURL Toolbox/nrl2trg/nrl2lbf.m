function [bsrf, du]=nrl2lbf(srf)

% Transform a nurls triangle into a Lagrange blending funaction triangle
% 
% Calling Sequences:
%
%     [bsrf, du]=nrl2lbf(srf)
%
% INPUTS:
%
%      srf - a nurls surface
%
% OUTPUT:
% 
%     bsrf    :   a blending function triangle. If du<0, nothing is done
%                   and an empty variable is returned
%     du :   an index of the shape  of the patch or the direction of u axis
%               du = -1 means the patch is an ellipse 
%               du = 0 means the patch is a quadrangle 
%               du = 1 or 2 is the axis NOT pointed to the focus point 
% 

% Transform the surface to make v axis point to focus point
[srf, du]=trgsrfdirect(srf); 
if du<0
    bsrf=[];
    return;
end

% Get the coordinates of the first two corner points of the triangle
pts=[0, 1; 0, 0];
pnts = nrleval(srf, pts);
P1=pnts(:, 1); P2=pnts(:, 2); 

% Evaluate the blending function triangle surface
crvt = nrlextract(srf);
crvs=crvt([3, 2, 1]);
t1=srf.knots{1}; 
t2=srf.knots{2}; 
[u, v]=meshgrid(t1, t2);
S=u.*(1-v); T=v;  
ct=trianglcoef(S(:)', T(:)');
p1=nrleval(crvs(1), S(:)');
p2=nrleval(crvs(2), T(:)');
p3=nrleval(crvs(3), T(:)');
x=ct{1}.*p1(1,:)+ct{2}.*p2(1,:)+ct{3}.*p3(1,:)+ct{4}*P1(1)+ct{5}*P2(1);
y=ct{1}.*p1(2,:)+ct{2}.*p2(2,:)+ct{3}.*p3(2,:)+ct{4}*P1(2)+ct{5}*P2(2);
z=ct{1}.*p1(3,:)+ct{2}.*p2(3,:)+ct{3}.*p3(3,:)+ct{4}*P1(3)+ct{5}*P2(3);
pb=[x; y; z];
pnts = nrleval(srf, {t1, t2});
pnts=permute(pnts, [1 3 2]);
pn=reshape(pnts, 3, []);

% Get nuls approximation of error of the triangle (ep)
m=length(t1); n=length(t2);
ep=pn-pb; 
coefs=permute(reshape(ep, 3, n, m), [1 3 2]);
esrf=srf; 
esrf.coefs(1:3,:)=coefs(:,:);

% Make a blending function triangle
bsrf.form='L-LBF';
bsrf.faces=srf;
bsrf.edges=crvs;
bsrf.corners=[P1, P2];
bsrf.efaces=esrf;
bsrf.direct=du;




