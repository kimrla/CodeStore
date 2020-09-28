function [pnt, jac, geom, hess] = nuhsrfgintvdeval (srf, order, knots, tt)

% nulsrfgintvdeval: Evaluation weights of the first derivative of NURL plane surface in Cartesion coordinates.
% 
% Calling Sequences:
%
%     [pnt, jac, geom, hess] = nuhsrfgintvdeval (order, coefs, knots, tt)
%
% INPUTS:
%
%     srf - a nurl surface
%
%      order - order of the nurl basis 
% 
%      knots - knot vectors of an interval
% 
%      tt     - parametric evaluation points. 
%
% OUTPUT:
%
%   pnt  - interpolation matrix used to evaluate points.
%
%   jac  - interpolation matrix used to evaluate first derivatives (Jacobian).
%
%   geom - the values, derivatives and determine of Jacobian at the evaluation nodes
%
%   hess - interpolation matrix used to evaluate second derivatives (Hessian).
% 

dd=max(srf.coefs(:))-min(srf.coefs(:));
if max(srf.coefs(3,:))-min(srf.coefs(3,:))>(1e-6)*dd
    error('Input is not a plane surface.');
end

if (iscell(knots) && numel(knots)==2)
    num1 = 2*length(knots{1});
    num2 = 2*length(knots{2});

    % Get derivatives and weights of derivatives in parametric coordinates
    if nargout<=3
        [pt, jac] = nrldeval (srf, tt);
        dpu=jac{1}; dpv=jac{2}; 
        [pnt, jaci] = nuhgintvdeval (order, knots, tt);
        Jacb=dpu(1,:).*dpv(2,:)-dpu(2,:).*dpv(1,:);
        geom={Jacb, pt, jac};
    elseif nargout==4
        [pt, jac, hess] = nrldeval (srf, tt);
        dpu=jac{1}; dpv=jac{2}; 
        dpuu=hess{1,1}; dpuv=hess{1,2}; dpvv=hess{2,2}; 
        [pnt, jaci, hessi] = nuhgintvdeval (order, knots, tt);
        Jacb=dpu(1,:).*dpv(2,:)-dpu(2,:).*dpv(1,:);
        geom={Jacb, pt, jac, hess};
    end

    % Get first derivatives and weights of derivatives in Cartesian coordinates    
    dpxi=(jaci{1}.*repmat(dpv(2, :), [num1*num2, 1]) ...
            - jaci{2}.*repmat(dpu(2, :), [num1*num2, 1]))./ ...
              repmat(Jacb, [num1*num2, 1]); 
    dpyi=(jaci{2}.*repmat(dpu(1, :), [num1*num2, 1]) ...
            - jaci{1}.*repmat(dpv(1, :), [num1*num2, 1]))./ ...
              repmat(Jacb, [num1*num2, 1]);      
     jac={dpxi, dpyi};
     
     if nargout==4
         dpxxi=(repmat(dpv(2, :), [num1*num2, 1]).^2.*( hessi{1,1} ...
                 - repmat(dpuu(1, :), [num1*num2, 1]).*dpxi ...
                 - repmat(dpuu(2, :), [num1*num2, 1]).*dpyi ) ...
            + repmat(dpu(2, :), [num1*num2, 1]).^2.*(hessi{2,2} ...
                 - repmat(dpvv(1, :), [num1*num2, 1]).*dpxi ...
                 - repmat(dpvv(2, :), [num1*num2, 1]).*dpyi) ...
             -2*repmat(dpv(2, :), [num1*num2, 1]).*repmat(dpu(2, :), [num1*num2, 1]).*(hessi{1,2} ...
                 - repmat(dpuv(1, :), [num1*num2, 1]).*dpxi ...
                 - repmat(dpuv(2, :), [num1*num2, 1]).*dpyi))./ ...
              repmat(Jacb, [num1*num2, 1]).^2; 
          dpyyi=(repmat(dpv(1, :), [num1*num2, 1]).^2.*( hessi{1,1} ...
                 - repmat(dpuu(1, :), [num1*num2, 1]).*dpxi ...
                 - repmat(dpuu(2, :), [num1*num2, 1]).*dpyi ) ...
            + repmat(dpu(1, :), [num1*num2, 1]).^2.*(hessi{2,2} ...
                 - repmat(dpvv(1, :), [num1*num2, 1]).*dpxi ...
                 - repmat(dpvv(2, :), [num1*num2, 1]).*dpyi) ...
             -2*repmat(dpv(1, :), [num1*num2, 1]).*repmat(dpu(1, :), [num1*num2, 1]).*(hessi{1,2} ...
                 - repmat(dpuv(1, :), [num1*num2, 1]).*dpxi ...
                 - repmat(dpuv(2, :), [num1*num2, 1]).*dpyi))./ ...
              repmat(Jacb, [num1*num2, 1]).^2; 
          dpxyi=(-repmat(dpv(1, :), [num1*num2, 1]).*repmat(dpv(2, :), [num1*num2, 1]).*( hessi{1,1} ...
                 - repmat(dpuu(1, :), [num1*num2, 1]).*dpxi ...
                 - repmat(dpuu(2, :), [num1*num2, 1]).*dpyi ) ...
            - repmat(dpu(1, :), [num1*num2, 1]).*repmat(dpu(2, :), [num1*num2, 1]).*(hessi{2,2} ...
                 - repmat(dpvv(1, :), [num1*num2, 1]).*dpxi ...
                 - repmat(dpvv(2, :), [num1*num2, 1]).*dpyi) ...
             + (repmat(dpu(1, :), [num1*num2, 1]).*repmat(dpv(2, :), [num1*num2, 1]) ...
                  + repmat(dpu(2, :), [num1*num2, 1]).*repmat(dpv(1, :), [num1*num2, 1]) ).* ...
                  (hessi{1,2} - repmat(dpuv(1, :), [num1*num2, 1]).*dpxi ...
                 - repmat(dpuv(2, :), [num1*num2, 1]).*dpyi))./ ...
              repmat(Jacb, [num1*num2, 1]).^2; 
          hess={dpxxi, dpxyi; dpxyi, dpyyi};
     end
else
    error('NURL structure is not a surface.');
end

%% demo
% crv1 = nrlcirc(0.5, [0, 0], 0, pi/2);
% crv2 = nrlcirc(1, [0, 0], 0, pi/2);
% srf=nrlruled(crv1, crv2);
% srf = nrldegelev(srf, [2, 2]);
% srf=nrlkntins(srf, {15, 8});
% nrlplot(srf, [100, 100], 'ctrl');
% view(2); axis equal;
% 
% % Get first derivatives and weights of derivatives in Cartesian coordinates
% m=10; n=13;
% s=linspace(0, 1, m); t=linspace(0, 1, n);
% [intpi, jaci, ~, hessi] = nuhsrfgintvdeval (srf, srf.order, srf.knots, {s, t});
% dpxi=jaci{1}; dpyi=jaci{2}; 
% dpxxi=hessi{1,1}; 
% dpyyi=hessi{2,2}; 
% dpxyi=hessi{1,2}; 
% 
% % Prepare the coefficients vector
% fun = @ (x, y) sin(x).*sin(y);
% funx = @ (x, y) cos(x).*sin(y);
% funy = @ (x, y) sin(x).*cos(y);
% funxx = @ (x, y) -sin(x).*sin(y);
% funyy = @ (x, y) -sin(x).*sin(y);
% funxy = @ (x, y) cos(x).*cos(y);
% 
% % The points and values on knots
% [pnt0, jac, hess] = nrldeval (srf, srf.knots);
% x=squeeze(pnt0(1,:,:)); y=squeeze(pnt0(2,:,:)); 
% xs=squeeze(jac{1}(1,:,:)); ys=squeeze(jac{1}(2,:,:)); 
% xt=squeeze(jac{2}(1,:,:)); yt=squeeze(jac{2}(2,:,:)); 
% xst=squeeze(hess{1,2}(1,:,:)); yst=squeeze(hess{1,2}(2,:,:)); 
% F0=fun(x, y); Fx0=funx(x, y); 
% Fy0=funy(x, y); Fxx0=funxx(x, y); 
% Fyy0=funyy(x, y); Fxy0=funxy(x, y); 
% Fs0=Fx0.*xs+Fy0.*ys;
% Ft0=Fx0.*xt+Fy0.*yt;
% Fst0=Fxx0.*xs.*xt+Fxy0.*(xs.*yt+xt.*ys)+Fyy0.*ys.*yt+Fx0.*xst+Fy0.*yst;
% 
% % Get grid indexes
% num1 = length(srf.knots{1});
% num2 = length(srf.knots{2});
% [pt, jc, hs] = nuhindex(num1, num2);
% ht=zeros(1, 4*num1*num2);
% ht(pt)=F0(:);
% ht(jc{1})=Fs0(:);
% ht(jc{2})=Ft0(:);
% ht(hs)=Fst0(:);
% 
% % Test interpolation in Cartesian coordinates
% pnt = nrldeval (srf, {s, t});
% xi=squeeze(pnt(1,:,:)); 
% yi=squeeze(pnt(2,:,:)); 
% F=sin(xi).*sin(yi);
% Fx=cos(xi).*sin(yi);
% Fy=sin(xi).*cos(yi);
% Fxx=-sin(xi).*sin(yi);
% Fyy=-sin(xi).*sin(yi);
% Fxy=cos(xi).*cos(yi);
% 
% Fi=reshape(ht*intpi, m, n);
% figure; surf(xi, yi, Fi); title('Fi');
% figure; surf(xi, yi, Fi-F); title('error of Fi');
% 
% Fxi=reshape(ht*dpxi, m, n);
% figure; surf(xi, yi, Fxi); title('Fxi');
% figure; surf(xi, yi, Fxi-Fx); title('error of Fxi');
% 
% Fyi=reshape(ht*dpyi, m, n);
% figure; surf(xi, yi, Fyi); title('Fyi');
% figure; surf(xi, yi, Fyi-Fy); title('error of Fyi');
% 
% Fxxi=reshape(ht*dpxxi, m, n);
% figure; surf(xi, yi, Fxxi); title('Fxxi');
% figure; surf(xi, yi, Fxxi-Fxx); title('error of Fxxi');
% 
% Fyyi=reshape(ht*dpyyi, m, n);
% figure; surf(xi, yi, Fyyi); title('Fyyi');
% figure; surf(xi, yi, Fyyi-Fyy); title('error of Fyyi');
% 
% Fxyi=reshape(ht*dpxyi, m, n);
% figure; surf(xi, yi, Fxyi); title('Fxyi');
% figure; surf(xi, yi, Fxyi-Fxy); title('error of Fxyi');




