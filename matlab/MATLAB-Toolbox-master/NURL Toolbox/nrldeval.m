function varargout = nrldeval (nurl, tt)

% NRLDEVAL: Evaluation of the derivative and second derivatives of NURL curve, surface or volume.
% 
% Calling Sequences:
%
%     [pnt, jac] = nrldeval (crv, tt)
%     [pnt, jac] = nrldeval (srf, {tu tv})
%     [pnt, jac] = nrldeval (vol, {tu tv tw})
%     [pnt, jac, hess] = nrldeval (crv, tt)
%     [pnt, jac, hess] = nrldeval (srf, {tu tv})
%     [pnt, jac, hess] = nrldeval (vol, {tu tv tw})
%
% INPUTS:
%
%   crv,   srf,   vol   - original NURL curve, surface or volume.
%   tt     - parametric evaluation points
%            If the nurl is a surface or a volume then tt is a cell
%            {tu, tv} or {tu, tv, tw} are the parametric coordinates
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

[cp, cw] = nrleval(nurl, tt);

if (iscell(nurl.knots))
  if (size(nurl.knots,2) == 3)
  % NURL structure represents a volume
    temp = cw(ones(3,1),:,:,:);
    pnt = cp./temp;
  
    [cup,cuw] = nrleval (nurl, tt, [1, 0, 0]);
    tempu = cuw(ones(3,1),:,:,:);
    jac{1} = (cup-tempu.*pnt)./temp;
  
    [cvp,cvw] = nrleval (nurl, tt, [0, 1, 0]);
    tempv = cvw(ones(3,1),:,:,:);
    jac{2} = (cvp-tempv.*pnt)./temp;

    [cwp,cww] = nrleval (nurl, tt, [0, 0, 1]);
    tempw = cww(ones(3,1),:,:,:);
    jac{3} = (cwp-tempw.*pnt)./temp;

% second derivatives
    if (nargout == 3)
        [cuup, cuuw] = nrleval (nurl, tt, [2, 0, 0]);
        tempuu = cuuw(ones(3,1),:,:,:);
        hess{1,1} = (cuup - (2*cup.*tempu + cp.*tempuu)./temp + 2*cp.*tempu.^2./temp.^2)./temp;
        clear cuup cuuw tempuu

        [cvvp, cvvw] = nrleval (nurl, tt, [0, 2, 0]);
        tempvv = cvvw(ones(3,1),:,:,:);
        hess{2,2} = (cvvp - (2*cvp.*tempv + cp.*tempvv)./temp + 2*cp.*tempv.^2./temp.^2)./temp;
        clear cvvp cvvw tempvv

        [cwwp, cwww] = nrleval (nurl, tt, [0, 0, 2]);
        tempww = cwww(ones(3,1),:,:,:);
        hess{3,3} = (cwwp - (2*cwp.*tempw + cp.*tempww)./temp + 2*cp.*tempw.^2./temp.^2)./temp;
        clear cwwp cwww tempww

        [cuvp, cuvw] = nrleval (nurl, tt, [1, 1, 0]);
        tempuv = cuvw(ones(3,1),:,:,:);
        hess{1,2} = (cuvp - (cup.*tempv + cvp.*tempu + cp.*tempuv)./temp + 2*cp.*tempu.*tempv./temp.^2)./temp;
        hess{2,1} = hess{1,2};
        clear cuvp cuvw tempuv

        [cuwp, cuww] = nrleval (nurl, tt, [1, 0, 1]);
        tempuw = cuww(ones(3,1),:,:,:);
        hess{1,3} = (cuwp - (cup.*tempw + cwp.*tempu + cp.*tempuw)./temp + 2*cp.*tempu.*tempw./temp.^2)./temp;
        hess{3,1} = hess{1,3};
        clear cuwp cuww tempuw

        [cvwp, cvww] = nrleval (nurl, tt, [0, 1, 1]);
        tempvw = cvww(ones(3,1),:,:,:);
        hess{2,3} = (cvwp - (cvp.*tempw + cwp.*tempv + cp.*tempvw)./temp + 2*cp.*tempv.*tempw./temp.^2)./temp;
        hess{3,2} = hess{2,3};
        clear cvwp cvww tempvw
    end

  elseif (size(nurl.knots,2) == 2)
    % NURL structure represents a surface 
    temp = cw(ones(3,1),:,:);
    pnt = cp./temp;
  
    [cup, cuw] = nrleval (nurl, tt, [1, 0]);
    tempu = cuw(ones(3,1),:,:);
    jac{1} = (cup-tempu.*pnt)./temp;
  
    [cvp, cvw] = nrleval (nurl, tt, [0, 1]);
    tempv = cvw(ones(3, 1), :, :);
    jac{2} = (cvp-tempv.*pnt)./temp;

% second derivatives
    if (nargout == 3) 
        [cuup, cuuw] = nrleval (nurl, tt, [2, 0]);
        tempuu = cuuw(ones(3,1),:,:);
        hess{1,1} = (cuup - (2*cup.*tempu + cp.*tempuu)./temp + 2*cp.*tempu.^2./temp.^2)./temp;

        [cvvp, cvvw] = nrleval (nurl, tt, [0, 2]);
        tempvv = cvvw(ones(3,1),:,:);
        hess{2,2} = (cvvp - (2*cvp.*tempv + cp.*tempvv)./temp + 2*cp.*tempv.^2./temp.^2)./temp;

        [cuvp, cuvw] = nrleval (nurl, tt, [1, 1]);
        tempuv = cuvw(ones(3,1),:,:);
        hess{1,2} = (cuvp - (cup.*tempv + cvp.*tempu + cp.*tempuv)./temp + 2*cp.*tempu.*tempv./temp.^2)./temp;
        hess{2,1} = hess{1,2};
    end
  end
else
  % NURL is a curve  
  temp = cw(ones(3,1),:);
  pnt = cp./temp;
  
  % first derivative
  [cup, cuw]=nrleval (nurl, tt, 1);
  temp1 = cuw(ones(3,1),:);
  jac = (cup-temp1.*pnt)./temp;
  if (iscell (tt))
    jac = {jac};
  end

  % second derivative
  if (nargout == 3)
    [cuup, cuuw] = nrleval (nurl, tt, 2);
    temp2 = cuuw(ones(3,1),:);
    hess = (cuup - (2*cup.*temp1 + cp.*temp2)./temp + 2*cp.*temp1.^2./temp.^2)./temp;
    if (iscell (tt))
      hess = {hess};
    end
  end
  
end

varargout{1} = pnt;
varargout{2} = jac;
if (nargout == 3)
  varargout{3} = hess;
end

end

%% Demo - curve
% a=2; b=1; N=6;
% sang=0; eang=2*pi;
% center=[0, 0];
% 
% crv = nrlellip(a, b, center, sang, eang);
% 
% s=linspace(0, 1, 200); 
% ang=(eang-sang)*s+sang;
% x=a*cos(ang)+center(1);
% y=b*sin(ang)+center(2);
% 
% t=linspace(0, 1, 13); 
% [p1, dp1, ddp] = nrldeval(crv, t);
% figure; hold on;
% plot(x, y); 
% plot(p1(1,:),p1(2,:),'ro');
% h = quiver(p1(1,:),p1(2,:),dp1(1,:),dp1(2,:));
% axis equal;

%% Demo - surface
% R=1; N=6;
% s1=0; s2=pi; t1=pi/3; t2=2*pi/3;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% figure; hold on;
% nrlplot(srf, [100, 100], 'ctrl');
% nrlplot(srf, [10, 10], 'quiver');
% view(3); axis equal;

%% Demo - Volume
% R=1; N=6;
% s1=0; s2=pi/2; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;
% 
% vol  = nrlextrude (srf, [0.2 0.2 0.2]);
% 
% figure; hold on;
% nrlplot(vol, [20, 21, 8], 'ctrl');
% nrlplot(vol, [5, 6, 3], 'quiver');
% view(3); axis equal;



