function [pnt, jac, geom, hess] = nrlsrfgintvdeval (order, coefs, knots, tt)

% nrlsrfgintvdeval: Evaluation weights of the first derivative of NURL plane surface in Cartesion coordinates.
% 
% Calling Sequences:
%
%     [pnt, jac, geom, hess] = nrlsrfgintvdeval (order, coefs, knots, tt)
%
% INPUTS:
%
%      order - order of the nurl basis 
% 
%      coefs -  Cartesian coordinates and weights
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

dd=max(coefs(:))-min(coefs(:));
if max(coefs(3,:))-min(coefs(3,:))>(1e-6)*dd
    error('Input is not a plane surface.');
end

if (iscell(knots) && numel(knots)==2)
    num1 = length(knots{1});
    num2 = length(knots{2});

    % Get derivatives and weights of derivatives in parametric coordinates
    weights=squeeze(coefs(4,:,:));
    coefs=reshape(coefs(1:3,:), 3, []);
    if nargout<=3
        [pnt, jaci] = nrlgintvdeval (order, weights, knots, tt);
        pt=coefs*pnt;
        dpu=coefs*jaci{1}; 
        dpv=coefs*jaci{2}; 
        Jacb=dpu(1,:).*dpv(2,:)-dpu(2,:).*dpv(1,:);
        geom={Jacb, pt, {dpu, dpv}};
    elseif nargout==4
        [pnt, jaci, hessi] = nrlgintvdeval (order, weights, knots, tt);
        pt=coefs*pnt;
        dpu=coefs*jaci{1}; 
        dpv=coefs*jaci{2}; 
        dpuu=coefs*hessi{1,1}; 
        dpuv=coefs*hessi{1,2}; 
        dpvv=coefs*hessi{2,2}; 
        Jacb=dpu(1,:).*dpv(2,:)-dpu(2,:).*dpv(1,:);
        geom={Jacb, pt, {dpu, dpv}, {dpuu, dpuv; dpuv, dpvv}};
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


%% Demo
% crv1 = nrlcirc(0.5, [0, 0], 0, pi/2);
% crv2 = nrlcirc(1, [0, 0], 0, pi/2);
% srf=nrlruled(crv1, crv2);
% srf = nrldegelev(srf, [2, 2]);
% srf=nrlkntins(srf, {25, 15});
% nrlplot(srf, [100, 100], 'ctrl');
% view(2); axis equal;
% 
% % Get first derivatives and weights of derivatives in Cartesian coordinates
% m=10; n=13;
% s=linspace(0, 1, m); t=linspace(0, 1, n);
% [intpi, jac, ~, hess] = nrlsrfgintvdeval (srf.order, srf.coefs, srf.knots, {s, t});
% dpxi=jac{1}; dpyi=jac{2}; 
% dpxxi=hess{1,1}; 
% dpyyi=hess{2,2}; 
% dpxyi=hess{1,2}; 
% 
% % Test interpolation in Cartesian coordinates
% num1 = length(srf.knots{1});
% num2 = length(srf.knots{2});
% x=reshape(srf.coefs(1, :), num1, num2);
% y=reshape(srf.coefs(2, :), num1, num2);
% F0=sin(x).*sin(y);
% 
% xi=reshape(x(:)'*intpi, m, n);
% yi=reshape(y(:)'*intpi, m, n);
% F=sin(xi).*sin(yi);
% Fx=cos(xi).*sin(yi);
% Fy=sin(xi).*cos(yi);
% Fxx=-sin(xi).*sin(yi);
% Fyy=-sin(xi).*sin(yi);
% Fxy=cos(xi).*cos(yi);
% 
% Fi=reshape(F0(:)'*intpi, m, n);
% figure; surf(xi, yi, Fi); title('Fi');
% figure; surf(xi, yi, Fi-F); title('error of Fi');
% 
% Fxi=reshape(F0(:)'*dpxi, m, n);
% figure; surf(xi, yi, Fxi); title('Fxi');
% figure; surf(xi, yi, Fxi-Fx); title('error of Fxi');
% 
% Fyi=reshape(F0(:)'*dpyi, m, n);
% figure; surf(xi, yi, Fyi); title('Fyi');
% figure; surf(xi, yi, Fyi-Fy); title('error of Fyi');
% 
% Fxxi=reshape(F0(:)'*dpxxi, m, n);
% figure; surf(xi, yi, Fxxi); title('Fxxi');
% figure; surf(xi, yi, Fxxi-Fxx); title('error of Fxxi');
% 
% Fyyi=reshape(F0(:)'*dpyyi, m, n);
% figure; surf(xi, yi, Fyyi); title('Fyyi');
% figure; surf(xi, yi, Fyyi-Fyy); title('error of Fyyi');
% 
% Fxyi=reshape(F0(:)'*dpxyi, m, n);
% figure; surf(xi, yi, Fxyi); title('Fxyi');
% figure; surf(xi, yi, Fxyi-Fxy); title('error of Fxyi');








