function varargout = nrlgintvdeval (order, weights, knots, tt)

% nrlgintvdeval: Evaluation matrices of first and second derivatives of NURL curve, surface or volume.
% 
% Calling Sequences:
%
%     [pnti, jaci] = nrlgintvdeval (order, weights, knots, tt)
%     [pnti, jaci] = nrlgintvdeval (order, weights, knots, {tu tv})
%     [pnti, jaci] = nrlgintvdeval (order, weights, knots, {tu tv tw})
%     [pnti, jaci, hessi] = nrlgintvdeval (order, weights, knots, tt)
%     [pnti, jaci, hessi] = nrlgintvdeval (order, weights, knots, {tu tv})
%     [pnti, jaci, hessi] = nrlgintvdeval (order, weights, knots, {tu tv tw})
%
% INPUTS:
%
%      order - order of the nurls basis 
% 
%      weights -  Cartesian weights of an interval
% 
%      knots - knot vectors of an interval
% 
%      tt     - parametric evaluation points. If the nurls is 
%            a surface or a volume then tt is a cell
%            {tu, tv} or {tu, tv, tw} are the parametric coordinates
%
% OUTPUT:
%
%   pnti  - interpolation matrix used to evaluate points.
%   jaci  - interpolation matrix used to evaluate first derivatives (Jacobian).
%   hessi - interpolation matrix used to evaluate second derivatives (Hessian).
%
% See also:
%  
%     nrlgdeval
%

[cp, cw] = nrlgintveval(order, weights, knots, tt);

if (iscell(knots))
  if (size(knots,2) == 3)
  % NURL structure represents a volume
    num1 = length(knots{1});
    num2 = length(knots{2});
    num3 = length(knots{3});
    temp = repmat(reshape(cw, 1, []), [num1*num2*num3, 1]); 
    pnti = cp./temp;
  
    [cup,cuw] = nrlgintveval(order, weights, knots, tt, [1, 0, 0]);
    tempu = repmat(reshape(cuw, 1, []), [num1*num2*num3, 1]); 
    jaci{1} = (cup-tempu.*pnti)./temp;
  
    [cvp,cvw] = nrlgintveval(order, weights, knots, tt, [0, 1, 0]);
    tempv = repmat(reshape(cvw, 1, []), [num1*num2*num3, 1]); 
    jaci{2} = (cvp-tempv.*pnti)./temp;

    [cwp,cww] = nrlgintveval(order, weights, knots, tt, [0, 0, 1]);
    tempw = repmat(reshape(cww, 1, []), [num1*num2*num3, 1]); 
    jaci{3} = (cwp-tempw.*pnti)./temp;

% second derivatives
    if (nargout == 3)
        [cuup, cuuw] = nrlgintveval(order, weights, knots, tt, [2, 0, 0]);
        tempuu = repmat(reshape(cuuw, 1, []), [num1*num2*num3, 1]); 
        hessi{1,1} = (cuup - (2*cup.*tempu + cp.*tempuu)./temp + 2*cp.*tempu.^2./temp.^2)./temp;
        clear cuup cuuw tempuu

        [cvvp, cvvw] = nrlgintveval(order, weights, knots, tt, [0, 2, 0]);
        tempvv = repmat(reshape(cvvw, 1, []), [num1*num2*num3, 1]); 
        hessi{2,2} = (cvvp - (2*cvp.*tempv + cp.*tempvv)./temp + 2*cp.*tempv.^2./temp.^2)./temp;
        clear cvvp cvvw tempvv

        [cwwp, cwww] = nrlgintveval(order, weights, knots, tt, [0, 0, 2]);
        tempww = repmat(reshape(cwww, 1, []), [num1*num2*num3, 1]); 
        hessi{3,3} = (cwwp - (2*cwp.*tempw + cp.*tempww)./temp + 2*cp.*tempw.^2./temp.^2)./temp;
        clear cwwp cwww tempww

        [cuvp, cuvw] = nrlgintveval(order, weights, knots, tt, [1, 1, 0]);
        tempuv = repmat(reshape(cuvw, 1, []), [num1*num2*num3, 1]); 
        hessi{1,2} = (cuvp - (cup.*tempv + cvp.*tempu + cp.*tempuv)./temp + 2*cp.*tempu.*tempv./temp.^2)./temp;
        hessi{2,1} = hessi{1,2};
        clear cuvp cuvw tempuv

        [cuwp, cuww] = nrlgintveval(order, weights, knots, tt, [1, 0, 1]);
        tempuw = repmat(reshape(cuww, 1, []), [num1*num2*num3, 1]); 
        hessi{1,3} = (cuwp - (cup.*tempw + cwp.*tempu + cp.*tempuw)./temp + 2*cp.*tempu.*tempw./temp.^2)./temp;
        hessi{3,1} = hessi{1,3};
        clear cuwp cuww tempuw

        [cvwp, cvww] = nrlgintveval(order, weights, knots, tt, [0, 1, 1]);
        tempvw = repmat(reshape(cvww, 1, []), [num1*num2*num3, 1]); 
        hessi{2,3} = (cvwp - (cvp.*tempw + cwp.*tempv + cp.*tempvw)./temp + 2*cp.*tempv.*tempw./temp.^2)./temp;
        hessi{3,2} = hessi{2,3};
        clear cvwp cvww tempvw
    end

  elseif (size(knots,2) == 2)
    % NURL structure represents a surface 
    num1 = length(knots{1});
    num2 = length(knots{2});
    temp = repmat(reshape(cw, 1, []), [num1*num2, 1]);
    pnti = cp./temp;
  
    [cup, cuw] = nrlgintveval(order, weights, knots, tt, [1, 0]);
    tempu = repmat(reshape(cuw, 1, []), [num1*num2, 1]);
    jaci{1} = (cup-tempu.*pnti)./temp;
  
    [cvp, cvw] = nrlgintveval(order, weights, knots, tt, [0, 1]);
    tempv = repmat(reshape(cvw, 1, []), [num1*num2, 1]);
    jaci{2} = (cvp-tempv.*pnti)./temp;

% second derivatives
    if (nargout == 3) 
        [cuup, cuuw] = nrlgintveval(order, weights, knots, tt, [2, 0]); 
        tempuu = repmat(reshape(cuuw, 1, []), [num1*num2, 1]); 
        hessi{1,1} = (cuup - (2*cup.*tempu + cp.*tempuu)./temp + 2*cp.*tempu.^2./temp.^2)./temp;

        [cvvp, cvvw] = nrlgintveval(order, weights, knots, tt, [0, 2]); 
        tempvv = repmat(reshape(cvvw, 1, []), [num1*num2, 1]); 
        hessi{2,2} = (cvvp - (2*cvp.*tempv + cp.*tempvv)./temp + 2*cp.*tempv.^2./temp.^2)./temp;

        [cuvp, cuvw] = nrlgintveval(order, weights, knots, tt, [1, 1]); 
        tempuv = repmat(reshape(cuvw, 1, []), [num1*num2, 1]); 
        hessi{1,2} = (cuvp - (cup.*tempv + cvp.*tempu + cp.*tempuv)./temp + 2*cp.*tempu.*tempv./temp.^2)./temp;
        hessi{2,1} = hessi{1,2};
    end
  end
else
  % NURL is a curve  
  temp = repmat(cw,[length(knots), 1]);
  pnti = cp./temp;
  
  % first derivative
  [cup, cuw]=nrlgintveval(order, weights, knots, tt, 1);
  temp1 = repmat(cuw,[length(knots), 1]);
  jaci = (cup-temp1.*pnti)./temp;
  
  if (iscell (tt))
    jaci = {jaci};
  end

  % second derivative
  if (nargout == 3)
    [cuup, cuuw] = nrlgintveval(order, weights, knots, tt, 2);
    temp2 = repmat(cuuw,[length(knots), 1]);
    hessi = (cuup - (2*cup.*temp1 + cp.*temp2)./temp + 2*cp.*temp1.^2./temp.^2)./temp;
    if (iscell (tt))
      hessi = {hessi};
    end
  end
  
end

varargout{1} = pnti;
varargout{2} = jaci;
if (nargout == 3)
  varargout{3} = hessi;
end

end


%% ! Demo - curve
% % Major and minor semi-axes (a, b)
% % Start and end angles (sang, eang)
% % Coordinates of the center (center)
% % The number of knots (N)
% a=2; b=1; N=6;
% sang=0; eang=pi;
% center=[0, 0];
% 
% % Get elliptic arcs
% crv = nrlellip(a, b, center, sang, eang);
% nrlplot(crv, 30); hold on;
% nrlplot(crv, 10, 'quiver');
% axis equal; 
% 
% % Test nrlgintveval
% n=10;
% t=linspace(0, 1, n);
% crv=nrlkntins(crv, 6);
% G = nrlgintveval(crv.order, crv.coefs(4,:), crv.knots, t);
% p=crv.coefs(1:3,:)*G;
% plot(p(1,:), p(2,:), 'ro');
% 
% % Test nrlgintvdeval
% [intp, jac, hess] = nrlgintvdeval(crv.order, crv.coefs(4,:), crv.knots, t);
% p1=crv.coefs(1:3,:)*intp;
% dp1=crv.coefs(1:3,:)*jac;
% ddp1=crv.coefs(1:3,:)*hess;
% figure; hold on;
% nrlplot(crv, 50); hold on;
% plot(p1(1,:), p1(2,:), 'ro');
% quiver(p1(1,:), p1(2,:), dp1(1,:), dp1(2,:));
% axis equal; 
% 
% n=50; t=linspace(0, 1, n);
% [intp, jac, hess] = nrlgintvdeval(crv.order, crv.coefs(4,:), crv.knots, t);
% dp2=crv.coefs(1:3,:)*jac;
% figure; hold on;
% plot(dp2(1,:), dp2(2,:));
% plot(dp1(1,:), dp1(2,:), 'ro');
% quiver(dp1(1,:), dp1(2,:), ddp1(1,:), ddp1(2,:));
% axis equal; 

%% demo - surface
% % Creat a sphere
% R=1; N=6;
% s1=0; s2=pi/2; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% srf=nrlkntins(srf, {2, 1});
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;
% 
% % Test nrlgintveval
% m=10; n=13;
% s=linspace(0, 1, m);
% t=linspace(0, 1, n);
% tt={s, t};
% knots=srf.knots;
% weights=squeeze(srf.coefs(4,:,:));
% der=[0, 0];
% order=srf.order;
% 
% num1 = length(knots{1});
% num2 = length(knots{2});
% degree = order;
% 
% nt1 = length(tt{1});
% nt2 = length(tt{2});
% 
% G = nrlgintveval(order, weights, knots, tt, der);
% coefs=reshape(srf.coefs(1:3,:), 3, []);
% p=coefs*G;
% hold on;
% plot3(p(1,:), p(2,:), p(3,:), 'ro');
% axis equal;
% 
% x=reshape(p(1,:), nt1, nt2);
% y=reshape(p(2,:), nt1, nt2);
% z=reshape(p(3,:), nt1, nt2);
% figure;
% nrlplot(srf, [30, 30], 'ctrl');
% hold on
% surf(x,y,z);
% axis equal;
% 
% % Test nrldgeval
% [G1, jac, hess] = nrlgintvdeval (order, weights, knots, tt);
% pt=coefs*G1; dpu=coefs*jac{1}; dpv=coefs*jac{2}; 
% figure;
% nrblplot(srf, [40 40], 'light', 'on');
% hold on;
% quiver3(pt(1,:), pt(2,:), pt(3,:), dpu(1,:), dpu(2,:), dpu(3,:))
% quiver3(pt(1,:), pt(2,:), pt(3,:), dpv(1,:), dpv(2,:), dpv(3,:))
% 
%  dpuu=coefs*hess{1,1}; dpuv=coefs*hess{1,2}; 
%  dpvv=coefs*hess{2,2}; 
% x=reshape(dpv(1,:), nt1, nt2);
% y=reshape(dpv(2,:), nt1, nt2);
% z=reshape(dpv(3,:), nt1, nt2);
% figure;
% surf(x,y,z);
% hold on
% quiver3(dpv(1,:), dpv(2,:), dpv(3,:), dpvv(1,:), dpvv(2,:), dpvv(3,:));
% quiver3(dpv(1,:), dpv(2,:), dpv(3,:), dpuv(1,:), dpuv(2,:), dpuv(3,:));
% axis equal;

%% demo - volume
% % Creat a sphere
% R=1; N=6;
% s1=0; s2=pi/2; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% vol  = nrlextrude (srf, [0.2 0.2 0.2]);
% vol=nrlkntins(vol, {2, 1, 1});
% figure;
% nrlplot(vol, [20, 21, 8], 'ctrl');
% view(3); axis equal;
% 
% % Test nrlgintveval
% m=10; n=13; k=4;
% t1=linspace(0, 1, m);
% t2=linspace(0, 1, n);
% t3=linspace(0, 1, k);
% tt={t1, t2, t3};
% knots=vol.knots;
% weights=squeeze(vol.coefs(4,:,:,:));
% der=[0, 0, 0];
% order=vol.order;
% 
% num1 = length(knots{1});
% num2 = length(knots{2});
% num3 = length(knots{3});
% degree = order;
% 
% nt1 = length(tt{1});
% nt2 = length(tt{2});
% nt3 = length(tt{3});
% 
% % Get and plot results
% G = nrlgintveval(order, weights, knots, tt, der);
% coefs=reshape(vol.coefs(1:3,:), 3, []); 
% p=coefs*G; 
% figure; hold on; 
% nrblplot(vol, [20 20 8], 'light', 'on'); 
% plot3(p(1,:), p(2,:), p(3,:), 'ro'); 
% view(3); axis equal; 
% 
% X=reshape(p(1,:), nt1, nt2, nt3);
% Y=reshape(p(2,:), nt1, nt2, nt3);
% Z=reshape(p(3,:), nt1, nt2, nt3);
% figure;
% nrlplot(vol, [20, 20, 6], 'ctrl');
% hold on
% m=nt1; n=nt2; k=nt3;
% surf(squeeze(X(1,:,:)), squeeze(Y(1,:,:)), squeeze(Z(1,:,:))); 
% surf(squeeze(X(m,:,:)), squeeze(Y(m,:,:)), squeeze(Z(m,:,:))); 
% surf(squeeze(X(:,1,:)), squeeze(Y(:,1,:)), squeeze(Z(:,1,:))); 
% surf(squeeze(X(:,n,:)), squeeze(Y(:,n,:)), squeeze(Z(:,n,:))); 
% surf(squeeze(X(:,:,1)), squeeze(Y(:,:,1)), squeeze(Z(:,:,1))); 
% surf(squeeze(X(:,:,k)), squeeze(Y(:,:,k)), squeeze(Z(:,:,k))); 
% axis equal;
% 
% [pnti, jaci, hessi] = nrlgintvdeval (order, weights, knots, tt);
% pt=coefs*pnti; 
% dpu=coefs*jaci{1}; dpv=coefs*jaci{2}; dpw=coefs*jaci{3}; 
% ddpuu=coefs*hessi{1,1}; ddpuv=coefs*hessi{1,2}; ddpuw=coefs*hessi{1,3}; 
% ddpvv=coefs*hessi{2,2}; ddpvw=coefs*hessi{2,3}; 
% ddpww=coefs*hessi{3,3};
% quiver3(pt(1,:), pt(2,:), pt(3,:), dpu(1,:), dpu(2,:), dpu(3,:));
% quiver3(pt(1,:), pt(2,:), pt(3,:), dpv(1,:), dpv(2,:), dpv(3,:));
% quiver3(pt(1,:), pt(2,:), pt(3,:), dpw(1,:), dpw(2,:), dpw(3,:));






