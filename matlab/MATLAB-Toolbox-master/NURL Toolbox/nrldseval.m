function dpnts = nrldseval (nurl, tt, der)

% NRLDEVAL: Evaluation of the derivatives of NURL curve, surface or volume.
%
%     dpnts = nrldseval (crv, tt, der)
%     dpnts = nrldseval (srf, {tu tv}, der)
%     dpnts = nrldseval (vol, {tu tv tw}, der)
%
% INPUTS:
%
%   crv,   srf,   vol   - original NURL curve, surface or volume.
%   tt     - parametric evaluation points
%            If the nurbs is a surface or a volume then tt is a cell
%            {tu, tv} or {tu, tv, tw} are the parametric coordinates
%   der - orders of derivatives to be evaluated for NURL. 
%
% OUTPUT:
%
%   pnt  - evaluated points.
%   jac  - evaluated first derivatives (Jacobian).
%   hess - evaluated second derivatives (Hessian).
%

if (nargin < 2)
  error ('nrlrdeval: wrong number of input parameters')
end

if (~isstruct(nurl))
  error('NURL representation is not structure!');
end

if (~strcmp(nurl.form,'L-NURL'))
  error('Not a recognised NURL representation');
end

% Define outputs
if iscell(nurl.knots)
  dpnts=cell(der+1); 
  dA=cell(der+1); 
  dw=cell(der+1); 
else
  dpnts=cell(1, der+1); 
  dA=cell(1, der); 
  dw=cell(1, der+1); 
end
  
if (iscell(nurl.knots))
  if (size(nurl.knots,2) == 3)
    % NURL structure represents a volume
    for p=1:der(3)+1
       for l=1:der(2)+1
            for k=1:der(1)+1
                [cp, cw]=nrleval (nurl, tt, [k-1, l-1, p-1]); 
                dw{k, l, p}=cw(ones(3,1),:,:,:);
                dA{k, l, p}=cp; 
                dpu=0; dpv=0; dpw=0; dpuvw=0; 
                for i=2:k 
                    dpu=dpu+nchoosek(k-1, i-1)*dw{i, 1, 1}.*dpnts{k-i+1, l, p}; 
                end
                for j=2:l 
                    dpv=dpv+nchoosek(l-1, j-1)*dw{1, j, 1}.*dpnts{k, l-j+1, p}; 
                end
                for m=2:p 
                    dpw=dpw+nchoosek(p-1, m-1)*dw{1, 1, m}.*dpnts{k, l, p-m+1}; 
                end
                for i=2:k 
                    for j=2:l 
                        for m=2:p
                            dpuvw=dpuvw+nchoosek(k-1, i-1)*nchoosek(l-1, j-1)*nchoosek(p-1, m-1)*dw{i, j, m}.*dpnts{k-i+1, l-j+1, p-m+1}; 
                        end
                    end
                end
                dpnts{k, l, p}=(dA{k, l, p}-dpu-dpv-dpw-dpuvw)./dw{1, 1, 1}; 
            end
       end
    end
  elseif (size(nurl.knots,2) == 2)
    % NURL structure represents a surface 
    for l=1:der(2)+1
        for k=1:der(1)+1
            [cp, cw]=nrleval (nurl, tt, [k-1, l-1]); 
            dw{k, l}=repmat(cw, 3, 1);
            dA{k, l}=cp; 
            dpu=0; dpv=0; dpuv=0; 
            for i=2:k 
                dpu=dpu+nchoosek(k-1, i-1)*dw{i, 1}.*dpnts{k-i+1, l}; 
            end
            for j=2:l 
                dpv=dpv+nchoosek(l-1, j-1)*dw{1, j}.*dpnts{k, l-j+1}; 
            end
            for i=2:k 
                for j=2:l 
                    dpuv=dpuv+nchoosek(k-1, i-1)*nchoosek(l-1, j-1)*dw{i, j}.*dpnts{k-i+1, l-j+1}; 
                end
            end
            dpnts{k, l}=(dA{k, l}-dpu-dpv-dpuv)./dw{1, 1}; 
        end
    end    
  end
else
  % NURL is a curve  
  [cp, cw] = nrleval(nurl, tt);
  dw{1}=cw(ones(3,1),:);
  dpnts{1}=cp./dw{1};
  
  % Derivatives
  for k=1:der
      [cup, cuw]=nrleval (nurl, tt, k); 
      dw{k+1}=cuw(ones(3, 1), :); 
      dA{k}=cup;
      dpu=0;
      for i=1:k
          dpu=dpu+nchoosek(k, i)*dw{i+1}.*dpnts{k-i+1};
      end
      dpnts{k+1}=(dA{k}-dpu)./dw{1};
  end
end

end

%% demo - curve
% a=2; b=1; N=6;
% sang=0; eang=2*pi;
% center=[0, 0];
% crv = nrlellip(a, b, center, sang, eang);
% 
% s=linspace(0, 1, 200); 
% ang=(eang-sang)*s+sang;
% x=a*cos(ang)+center(1);
% y=b*sin(ang)+center(2);
% 
% t=linspace(0, 1, 13); 
% dps = nrldseval(crv, t, 2);
% figure; hold on;
% plot(x, y); 
% plot(dps{1}(1,:), dps{1}(2,:),'ro');
% h = quiver(dps{1}(1,:),dps{1}(2,:),dps{2}(1,:),dps{2}(2,:));
% axis equal;

%% demo - surface
% R=1; N=6;
% s1=0; s2=pi; t1=pi/3; t2=2*pi/3;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;
% 
% m=20; n=25;
% t=linspace(0,1,m); k=linspace(0,1,n);
% dpnts = nrldseval (srf, {t, k}, [2,2]);
% figure;
% nrlplot(srf, [100, 100], 'ctrl');
% hold on;
% quiver3(dpnts{1,1}(1,:,:), dpnts{1,1}(2,:,:), dpnts{1,1}(3,:,:), ...
%     dpnts{1,2}(1,:,:), dpnts{1,2}(2,:,:), dpnts{1,2}(3,:,:)); 
% quiver3(dpnts{1,1}(1,:,:), dpnts{1,1}(2,:,:), dpnts{1,1}(3,:,:), ...
%     dpnts{2,1}(1,:,:), dpnts{2,1}(2,:,:), dpnts{2,1}(3,:,:)); 
% axis equal;

%% demo - surface
% R=1; N=6;
% s1=0; s2=pi/2; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;
% 
% m=5; n=6; k=3;
% ut=linspace(0, 1, m); 
% vt=linspace(0, 1, n); 
% wt=linspace(0, 1, k); 
% vol  = nrlextrude (srf, [0.2 0.2 0.2]);
% dpnts = nrldseval(vol, {ut, vt, wt}, [1, 1, 1]);    
% figure; 
% nrlplot(vol, [20, 21, 8], 'ctrl');
% hold on; 
% p=1; q=1; r=1;
% quiver3(dpnts{p,q,r}(1,:,:), dpnts{p,q,r}(2,:,:), dpnts{p,q,r}(3,:,:), ...
%     dpnts{2,1,1}(1,:,:), dpnts{2,1,1}(2,:,:), dpnts{2,1,1}(3,:,:)); 
% quiver3(dpnts{p,q,r}(1,:,:), dpnts{p,q,r}(2,:,:), dpnts{p,q,r}(3,:,:), ...
%     dpnts{1,2,1}(1,:,:), dpnts{1,2,1}(2,:,:), dpnts{1,2,1}(3,:,:)); 
% quiver3(dpnts{p,q,r}(1,:,:), dpnts{p,q,r}(2,:,:), dpnts{p,q,r}(3,:,:), ...
%     dpnts{1,1,2}(1,:,:), dpnts{1,1,2}(2,:,:), dpnts{1,1,2}(3,:,:)); 
% axis equal;




