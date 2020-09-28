function [G, Gs, Gt, S, T]=nrltrgmat(tsrf, tt)

% Get first derivatives of the nurl basis matrix in area coordinates of a triangle
% 
%  INPUT:
%
%    tsrf    -   a nurl triangle patch.
%
%      tt     - a cell {tu, tv} of parametric evaluation points. 
%                or scattered parametric coordinates (see nrleval).
%
% OUTPUT:
% 
%    G - a matrix of the basis
% 
%    Gs, Gt - matrices of first derivatives with respect to area coordinates
% 
%    S, T - area coordinates transformed from natural coodinates tt
%

if (nargin < 2)
  error('Not enough input arguments');
end

if (~isstruct(tsrf))
  error('NURL representation is not structure!');
end

if (~strcmp(tsrf.form,'T-NURL'))
  error('Not a recognised triangle NURL representation');
end

% Get basis on a unit quadrilatril and sum the end basis
srf=tsrf.faces;
num1=srf.number(1); num2=srf.number(2); 
TN=num1*(num2-1)+1;
order=srf.order; 
knots=srf.knots; 
weights=squeeze(srf.coefs(4,:,:));
[G, jac] = nrlgintvdeval(order, weights, knots, tt); 
Gu=jac{1}; Gv=jac{2};
G(TN,:)=sum(G(TN:end,:)); G(TN+1:end,:)=[];
Gu(TN,:)=sum(Gu(TN:end,:)); Gu(TN+1:end,:)=[];
Gv(TN,:)=sum(Gv(TN:end,:)); Gv(TN+1:end,:)=[];

% Cope with focus point
crvs=tsrf.edges; crv2=crvs(2); crv3=crvs(3);
[~, jac2] = nrlgintvdeval (crv2.order, crv2.coefs(4,:), crv2.knots, 1);
[~, jac3] = nrlgintvdeval (crv3.order, crv3.coefs(4,:), crv3.knots, 1);
gs=zeros(TN,1); gt=zeros(TN,1);
for j=1:num2
    p=(j-1)*num1+1;
    q=(j-1)*num1+num1;
    if j<num2
        gs(q)=-jac2(j);
        gs(p)=jac3(j);
        gt(p)=jac3(j);
    else
        gs(end)=jac3(j)-jac2(j);
        gt(p)=jac3(j);
    end
end

% Get the parametric points of a unit triangle
if iscell(tt)
    [v, u]=meshgrid(tt{2}, tt{1});
else
    [m,~]=size(tt);
    if m~=2
        tt=tt';
    end
    u=tt(1,:); v=tt(2,:);
end
S=u.*(1-v); T=v;  

% Get the basis on a unit triangle
sv=S(:)'; tv=T(:)';
pp=tv==1; qq=tv~=1; 
nr=length(find(pp));
gs=repmat(gs, 1, nr); gt=repmat(gt, 1, nr);
sv=repmat(sv, TN, 1);  tv=repmat(tv, TN, 1); 
pp=repmat(pp, TN, 1); qq=repmat(qq, TN, 1); 
Gs=zeros(size(Gu)); Gt=Gs;
Gs(pp)=gs(:,:); Gt(pp)=gt(:,:); 
Gs(qq)=Gu(qq)./(1-tv(qq));
Gt(qq)=sv(qq).*Gu(qq)./(1-tv(qq)).^2+Gv(qq);


%% Demo
% % The number integration nodes (m, n)
% m=10; n=12; 
% 
% % Creat a triangular nurls patch
% crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
% af=pi/3; a=cos(af); b=sin(af);
% crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
% crv3 = nrlcirc(1, [0, 0], 0, af);
% srf = nrltrgcoons(crv1, crv2, crv3);
% 
% % Transform the surface into  a nurls triangle patch
% srf=nrlkntins(srf, {2, 1});
% srf = nrldegelev(srf, [2, 2]);
% tsrf=nrl2trg(srf); 
% srf=tsrf.faces;
% 
% figure; hold on;
% nrlplot(srf, [20, 20]);
% nrlaxisplot(srf);
% view(2); axis equal;
% 
% % Coefficients of blending function on triangles
% v1 = srf.knots{1}; v2 = srf.knots{2}; 
% num1=length(v1); num2=length(v2); 
% NB=num1+2*num2-3; 
% TN=num1*(num2-1)+1;
% 
% % Displacement functions on triangles
% [t1, C1]=GaussLobattoR(m, 0, 1); 
% [t2, C2]=GaussLobattoR(n, 0, 1);  
% t1=linspace(0,1,m); t2=linspace(0,1,n); 
% tt={t1, t2}; 
% 
% % Basis on an triangle
% [G, Gs, Gt, S, T]=nrltrgmat(tsrf, tt);
% 
% [s1, r1]=meshgrid(v2(1:end-1), v1(1:end));
% pts=[r1(:)', 0; s1(:)', 1];
% pnts=nrleval(srf, pts);
% p=pnts*G;
% dps=pnts*Gs;
% dpt=pnts*Gt;
% 
% % Test for gridded points
% figure; hold on;
% x=reshape(p(1,:), m, n);
% y=reshape(p(2,:), m, n);
% z=reshape(p(3,:), m, n);
% surf(x,y,z);
% quiver(p(1,:), p(2,:), dps(1,:), dps(2,:));
% quiver(p(1,:), p(2,:), dpt(1,:), dpt(2,:));
% colormap summer;      
% axis equal;
% axis([-0.1, 1.1, -0.1,0.95]);
% 
% % Test for scattered points
% [u, v]=meshgrid(tt{1}, tt{2});
% [G, Gs, Gt, S, T]=nrltrgmat(tsrf, [u(:)'; v(:)']);
% p=pnts*G;
% dps=pnts*Gs;
% dpt=pnts*Gt;
% figure; hold on;
% surf(x,y,z);
% quiver(p(1,:), p(2,:), dps(1,:), dps(2,:));
% quiver(p(1,:), p(2,:), dpt(1,:), dpt(2,:));
% colormap summer;      
% axis equal;
% axis([-0.1, 1.1, -0.1,0.95]);





