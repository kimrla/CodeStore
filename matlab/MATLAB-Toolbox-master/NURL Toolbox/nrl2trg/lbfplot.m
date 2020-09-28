function lbfplot(bsrf, nuv, ctrl)

% Plot a lbf surface
% 
% Calling Sequences:
% 
%     lbfplot (srf, [nu, nv], ctrl)
%
%  Input: 
%
%    bsrf - a nurls surface
%    nuv - the number of points
%    ctrl - plot control points ('ctrl'), quiver ('quiver') 
%            or points ('pnts')
%

if nargin==2
    if ischar(nuv)
        ctrl=nuv;
        nuv=30*ones(size(bsrf.faces.number));
    else
        ctrl='non';
    end
elseif nargin==1
    ctrl='non';
    nuv=30*ones(size(nrls.number));    
end
tf = ishold;

if strcmp(ctrl, 'quiver')
    m=nuv(1); n=nuv(2);
    t1=linspace(0,1,m); t2=linspace(0,1,n); 
    tt={t1, t2};
    [pp, jac]=lbfdeval(bsrf, tt);
    dps=jac{1}; dpt=jac{2};

    x=reshape (pp(1,:), n, m);
    y=reshape (pp(2,:), n, m);
    z=reshape (pp(3,:), n, m);
    surf(x, y, z); hold on;
    quiver3(pp(1,:), pp(2,:), pp(3,:), dps(1,:), dps(2,:), dps(3,:));
    quiver3(pp(1,:), pp(2,:), pp(3,:), dpt(1,:), dpt(2,:), dpt(3,:));
    colormap summer;    
    title('Quiver plot of aera coordinates.');
    view(2); 
    axis equal;
    
    tt={t1, t2};
    [u, v]=meshgrid(tt{1}, tt{2});
    U=repmat(u(:)', 3, 1);
    V=repmat(v(:)', 3, 1);
    dpu=(1-V).*dps;
    dpv=-U.*dps+dpt;
    figure; 
    surf(x, y, z); hold on; 
    quiver3(pp(1,:), pp(2,:), pp(3,:), dpu(1,:), dpu(2,:), dpu(3,:));
    quiver3(pp(1,:), pp(2,:), pp(3,:), dpv(1,:), dpv(2,:), dpv(3,:));
    colormap summer; 
    title('Quiver plot of natural coordinates.');
    view(2);
elseif strcmp(ctrl, 'pnts')
    m=nuv(1); n=nuv(2);
    t1=linspace(0,1,m); t2=linspace(0,1,n);  
    pnts=lbfeval(bsrf, {t1, t2});
    plot3(pnts(1,:,:), pnts(2,:,:), pnts(3,:,:), 'ro');
else
    m=nuv(1); n=nuv(2);
    t1=linspace(0,1,m); t2=linspace(0,1,n);  
    pp=lbfeval(bsrf, {t1, t2});
    x=reshape (pp(1,:), n, m);
    y=reshape (pp(2,:), n, m);
    z=reshape (pp(3,:), n, m);
    surf(x, y, z);
    colormap summer;      
    if strcmp(ctrl, 'ctrl')
        shading interp;
        hold on;
        for ii = 1:size (bsrf.faces.coefs, 2)
            coefs = reshape (bsrf.faces.coefs(1:3,ii,:), 3, []);
            plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'k--')
            plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'r.','MarkerSize',20)
        end
        for jj = 1:size (bsrf.faces.coefs, 3)
            coefs = reshape (bsrf.faces.coefs(1:3,:,jj), 3, []);
            plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'k--');
            plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'r.','MarkerSize',20)
        end
    end
end
axis equal;

if tf
    hold on;
else
    hold off;
end



