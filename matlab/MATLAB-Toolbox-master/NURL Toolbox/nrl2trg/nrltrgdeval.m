function [pnt, jac]=nrltrgdeval(bsrf, tt)

% Evaluate the derivatives of an triganle with respect to area coodinates
% 
%  INPUT:
%
%    bsrf    :   a blending function triangle.
%
%      tt     - a cell {tu, tv} of parametric evaluation points. 
%                or scattered parametric coordinates (see nrleval).
%
% OUTPUT:
% 
%    pnt - points
% 
%    jac - first derivatives with respect to area coordinates
%

% Get derivatives of the quadrangle and the ends of two boundary curves
srf=bsrf.faces;
[pnt, jac] = nrldeval (srf, tt);
dpu=jac{1}; dpv=jac{2};
crvs=bsrf.edges; 
crv2=crvs(2); crv3=crvs(3);
[~, jac2] = nrldeval (crv2, 1);
[~, jac3] = nrldeval (crv3, 1);
ds=jac3-jac2;
dt=jac3;

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

% Get the derivatives on a unit triangle by axis transformation
if iscell(tt)
    m=length(tt{1}); n=length(tt{2});
    pp=T==1; qq=T~=1; 
    [p, q]=size(find(pp));
    ds=repmat(ds, 1, p, q);  
    dt=repmat(dt, 1, p, q);  
    sv=repmat(reshape(S, 1, m, n), 3, 1);  
    tv=repmat(reshape(T, 1, m, n), 3, 1); 
    pp=repmat(reshape(pp, 1, m, n), 3, 1);  
    qq=repmat(reshape(qq, 1, m, n), 3, 1); 
else
    nt=length(u); 
    pp=T==1; qq=T~=1; 
    p=length(find(pp));
    ds=repmat(ds, 1, p);  
    dt=repmat(dt, 1, p);  
    sv=repmat(reshape(S, 1, nt), 3, 1);  
    tv=repmat(reshape(T, 1, nt), 3, 1); 
    pp=repmat(reshape(pp, 1, nt), 3, 1);  
    qq=repmat(reshape(qq, 1, nt), 3, 1); 
end
dps=zeros(size(dpu)); dpt=dps;
dps(pp)=ds(:,:,:); dpt(pp)=dt(:,:,:); 
dps(qq)=dpu(qq)./(1-tv(qq));
dpt(qq)=sv(qq).*dpu(qq)./(1-tv(qq)).^2+dpv(qq);

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
% % Transform nurls surface into Lagrange blending surface
% bsrf=nrl2lbf(srf);
% 
% % Evaluate the derivatives of an triganle with respect to area coodinates
% t1=linspace(0,1,m); t2=linspace(0,1,n); 
% tt={t1, t2}; 
% [pnt, jac]=nrltrgdeval(bsrf, tt);
% dps=jac{1}; dpt=jac{2};
% 
% figure; hold on;
% x=squeeze(pnt(1,:,:));
% y=squeeze(pnt(2,:,:));
% z=squeeze(pnt(3,:,:));
% surf(x,y,z);
% quiver(pnt(1,:), pnt(2,:), dps(1,:), dps(2,:));
% quiver(pnt(1,:), pnt(2,:), dpt(1,:), dpt(2,:));
% colormap summer;      
% axis equal;
% axis([-0.1, 1.1, -0.1,0.95]);
% 
% [u, v]=meshgrid(tt{1}, tt{2});
% [pnt, jac]=nrltrgdeval(bsrf, [u(:)'; v(:)']);
% dps=jac{1}; dpt=jac{2};
% figure; hold on;
% surf(x,y,z);
% quiver(pnt(1,:), pnt(2,:), dps(1,:), dps(2,:));
% quiver(pnt(1,:), pnt(2,:), dpt(1,:), dpt(2,:));
% colormap summer;      
% axis equal;
% axis([-0.1, 1.1, -0.1,0.95]);







