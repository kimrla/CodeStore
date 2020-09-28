function nrlplot(nrl, nuv, ctrl)

% Plot a nurl curve, surface or volume (and their control points and quiver)
% 
% Calling Sequences:
% 
%     nrlplot (crv, nu, ctrl)
%     nrlplot (srf, [nu, nv], ctrl)
%     nrlplot (vol, [nu, nv, nw], ctrl)
%     nrlplot (crv, ctrl)
%     nrlplot (srf, ctrl)
%     nrlplot (vol, ctrl)
%
%  Input: 
% 
%    nrls - a nurl curve, surface or volume 
%    nuv - the number of points
%    ctrl - plot control points ('ctrl'), quiver ('quiver') 
%            or points ('pnts')
%

if nargin==2
    if ischar(nuv)
        ctrl=nuv;
        nuv=100*ones(size(nrl.number));
    else
        ctrl='non';
    end
elseif nargin==1
    ctrl='non';
    nuv=100*ones(size(nrl.number));    
end
tf = ishold;

trl=length(nrl.number);
if trl==3
    % NURL structure represents a volume
    m=nuv(1); n=nuv(2); k=nuv(3);
    ut=linspace(0, 1, m); 
    vt=linspace(0, 1, n); 
    wt=linspace(0, 1, k); 
    if strcmp(ctrl, 'quiver')
        [pnts, jac] = nrldeval(nrl, {ut, vt, wt});            
        quiver3(pnts(1,:,:), pnts(2,:,:), pnts(3,:,:), ...
            jac{1}(1,:,:), jac{1}(2,:,:), jac{1}(3,:,:)); 
        hold on;
        quiver3(pnts(1,:,:), pnts(2,:,:), pnts(3,:,:), ...
            jac{2}(1,:,:), jac{2}(2,:,:), jac{2}(3,:,:)); 
        quiver3(pnts(1,:,:), pnts(2,:,:), pnts(3,:,:), ...
            jac{3}(1,:,:), jac{3}(2,:,:), jac{3}(3,:,:)); 
    elseif strcmp(ctrl, 'pnts')
        pnts=nrleval(nrl, {ut, vt, wt});    
        plot3(pnts(1,:), pnts(2,:), pnts(3,:), 'ro');
    else
        pnts=nrleval(nrl, {ut, vt, wt});    
        X=zeros(m,n,k); Y=zeros(m,n,k); Z=zeros(m,n,k);
        for i=1:m
            for j=1:n
                for p=1:k
                    X(i,j,p)=pnts(1,i,j,p);
                    Y(i,j,p)=pnts(2,i,j,p);
                    Z(i,j,p)=pnts(3,i,j,p);
                end
            end
        end
        surf(squeeze(X(1,:,:)), squeeze(Y(1,:,:)), squeeze(Z(1,:,:))); hold on; 
        surf(squeeze(X(m,:,:)), squeeze(Y(m,:,:)), squeeze(Z(m,:,:))); 
        surf(squeeze(X(:,1,:)), squeeze(Y(:,1,:)), squeeze(Z(:,1,:))); 
        surf(squeeze(X(:,n,:)), squeeze(Y(:,n,:)), squeeze(Z(:,n,:))); 
        surf(squeeze(X(:,:,1)), squeeze(Y(:,:,1)), squeeze(Z(:,:,1))); 
        surf(squeeze(X(:,:,k)), squeeze(Y(:,:,k)), squeeze(Z(:,:,k))); 
        colormap summer;
        if strcmp(ctrl, 'ctrl')
            shading interp;
            hold on;
            plot3(nrl.coefs(1,:), nrl.coefs(2,:), nrl.coefs(3,:), 'r.', 'MarkerSize', 20); 
        end
    end
elseif trl==2
    % NURL structure represents a surface
    m=nuv(1); n=nuv(2);
    t=linspace(0,1,m); k=linspace(0,1,n);    
    if strcmp(ctrl, 'quiver')
        [pnts, jac] = nrldeval (nrl, {t, k});
        quiver3(pnts(1,:,:), pnts(2,:,:), pnts(3,:,:), ...
            jac{1}(1,:,:), jac{1}(2,:,:), jac{1}(3,:,:)); 
        hold on;
        quiver3(pnts(1,:,:), pnts(2,:,:), pnts(3,:,:), ...
            jac{2}(1,:,:), jac{2}(2,:,:), jac{2}(3,:,:)); 
    elseif strcmp(ctrl, 'pnts')
        pnts=nrleval(nrl, {t, k});    
        plot3(pnts(1,:), pnts(2,:), pnts(3,:), 'ro');
    else
        pnts=nrleval(nrl, {t, k});
        pntsx=zeros(m,n); pntsy=zeros(m,n); pntsz=zeros(m,n);
        for i=1:m
            for j=1:n
                pntsx(i,j)=pnts(1,i,j);
                pntsy(i,j)=pnts(2,i,j);
                pntsz(i,j)=pnts(3,i,j);
            end
        end
        surf(pntsx,pntsy,pntsz);
        colormap summer;        
        if strcmp(ctrl, 'ctrl')
            shading interp;
            hold on;
            plot3(nrl.coefs(1,:), nrl.coefs(2,:), nrl.coefs(3,:), 'r.', 'MarkerSize', 20); 
        end
    end
elseif trl==1
    % NURL structure represents a curve
    t=linspace(0, 1, nuv); 
    if strcmp(ctrl, 'quiver')
        [p1, dp1] = nrldeval(nrl, t);
        quiver3(p1(1,:),p1(2,:),p1(3,:),dp1(1,:),dp1(2,:),dp1(3,:));
        if max(p1(3,:))==min(p1(3,:))
            view(2);
        end
    elseif strcmp(ctrl, 'pnts')
        pnts=nrleval(nrl, t);    
        plot3(pnts(1,:), pnts(2,:), pnts(3,:), 'ro');
         if max(pnts(3,:))==min(pnts(3,:))
             view(2);
         end
    else
        pnts=nrleval(nrl, t); 
        plot3(pnts(1,:), pnts(2,:), pnts(3,:));
        if strcmp(ctrl, 'ctrl')
            plot3(nrl.coefs(1,:), nrl.coefs(2,:), nrl.coefs(3,:), 'r.', 'MarkerSize', 20); 
        end
        if max(pnts(3,:))==min(pnts(3,:))
             view(2);
        end
    end
end
axis equal;

if tf
    hold on;
else
    hold off;
end


