function [pnts, jac]=nrl2lbfdeval(bsrf, tt)

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
%     pnts    :   evaluated points corresponding to the 
%                 parametric points.
%     jac  : evaluated first derivatives (Jacobian).
% 


[pnts, jac] = nrldeval(bsrf.faces, tt);
pnts=permute(pnts, [1 3 2]);
pnts=reshape(pnts, 3, []);
dptu=permute(jac{1}, [1 3 2]);
dpu=reshape(dptu, 3, []);
dptv=permute(jac{2}, [1 3 2]);
dpv=reshape(dptv, 3, []);
du=bsrf.direct;
[S, T]=getst(tt, du);
ss = repmat(S(:)', 3, 1);  
tt = repmat(T(:)', 3, 1);
pp=tt==1; qq=tt~=1; 
dps=zeros(size(dpu)); dpt=dps; 
dpt(pp)=dpv(pp);
dps(qq)=dpu(qq)./(1-tt(qq));
dpt(qq)=ss(qq).*dpu(qq)./(1-tt(qq)).^2+dpv(qq);

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
% 
% % Quiver plot on area coordinates by transforming
% t1=linspace(0,1,m); t2=linspace(0,1,n); 
% tt={t1, t2};
% [pn, jac]=nrl2lbfdeval(bsrf, tt);
% dps=jac{1}; dpt=jac{2}; 
% 
% x=reshape (pn(1,:), n, m);
% y=reshape (pn(2,:), n, m);
% z=reshape (pn(3,:), n, m);
% figure;
% surf(x, y, z); hold on;
% quiver3(pn(1,:), pn(2,:), pn(3,:), dps(1,:), dps(2,:), dps(3,:));
% quiver3(pn(1,:), pn(2,:), pn(3,:), dpt(1,:), dpt(2,:), dpt(3,:));
% colormap summer;    
% title('Quiver plot of aera coordinates (transformed from natrual coordinates).');
% view(2); 
% axis equal;





