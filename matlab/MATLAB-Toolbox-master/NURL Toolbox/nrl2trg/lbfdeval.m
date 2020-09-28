function [pp, jac]=lbfdeval(bsrf, tt)

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
%                or scattered parametric coordinates (see nrleval).
%
% OUTPUT:
% 
%     pp    :   evaluated points corresponding to the 
%                 parametric points.
%     jac  : evaluated first derivatives (Jacobian).
% 

[u, v]=meshgrid(tt{1}, tt{2});
S=u.*(1-v); T=v; 
crvs=bsrf.edges;
P1=bsrf.corners(:,1);
P2=bsrf.corners(:,2);
ct=trianglcoef(S(:)', T(:)');
ctx=trianglcoefx(S(:)', T(:)');
cty=trianglcoefy(S(:)', T(:)');
[p1, dp1]=nrldeval(crvs(1), S(:)'); 
[p2, dp2]=nrldeval(crvs(2), T(:)'); 
[p3, dp3]=nrldeval(crvs(3), T(:)'); 
for i=1:5
    ct{i} = repmat(ct{i}, 3, 1);
    if i>3
       continue; 
    end
    ctx{i} = repmat(ctx{i}, 3, 1);
    cty{i} = repmat(cty{i}, 3, 1);
end
m=length(tt{1}); n=length(tt{2}); 
P1 = repmat(P1, 1, m*n); 
P2 = repmat(P2, 1, m*n);
ss = repmat(S(:)', 3, 1);  
pb=ct{1}.*p1+ct{2}.*p2+ct{3}.*p3+ct{4}.*P1+ct{5}.*P2;
dpbs=ct{1}.*dp1+ctx{1}.*(p2-p3)+ctx{2}.*(dp2 -dp3)+ctx{3}.*(P2-p1)+P1-P2;
dpbt=ct{2}.*dp2+ct{3}.*dp3+cty{1}.*(P2.*ss-p1)+cty{2}.*(P2-dp1)+cty{3}.*(p2-p3)+P1;

% Evaluate derivatives of the inners error surface
[pnts, jac] = nrldeval(bsrf.efaces, tt);
pnts=permute(pnts, [1 3 2]);
pa=reshape(pnts, 3, []);
dptu=permute(jac{1}, [1 3 2]);
dpu=reshape(dptu, 3, []);
dptv=permute(jac{2}, [1 3 2]);
dpv=reshape(dptv, 3, []);
tt = repmat(T(:)', 3, 1);
pp=tt==1; qq=tt~=1; 
dpas=zeros(size(dpbs)); dpat=dpas; 
dpat(pp)=dpv(pp);
dpas(qq)=dpu(qq)./(1-tt(qq));
dpat(qq)=ss(qq).*dpu(qq)./(1-tt(qq)).^2+dpv(qq);

% Assemble the results
pp=pb+pa;
dps=dpbs+dpas;
dpt=dpbt+dpat;
jac={dps, dpt};


%% ! Demo
% % The number integration nodes (m, n)
% m=10; n=12; 
% 
% % Creat a triangular nurls patch
% crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
% af=pi/3; a=cos(af); b=sin(af);
% crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
% crv3 = nrlcirc(1, [0, 0], 0, af);
% srf = nrltrgcoons(crv1, crv2, crv3);
% srf=nrlkntins(srf, {5, 4});
% srf = nrldegelev(srf, [2, 2]);
% 
% % Plot the nurls surface
% figure; hold on;
% nrlaxisplot(srf);
% view(2); axis equal;
% 
% % Transform nurls surface into Lagrange blending surface
% bsrf=nrl2lbf(srf);
% 
% % Evaluate points and first derivatives of the lbf surface
% figure;
% lbfplot(bsrf, [m, n], 'quiver');
% view(2);




