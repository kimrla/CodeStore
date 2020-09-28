function [pnt, jac, Jacb] = srftrgplaneder (tsrf, tt)

% srftrgplaneder: Evaluation weights of the first derivative of NURL triangular plane surface in Cartesian coordinates.
% 
% Calling Sequences:
%
%     [pnt, jac] = srftrgplaneder (order, coefs, knots, tt)
%
% INPUTS:
%
%    tsrf    -   a nurl triangle patch.
%
%      tt     - a cell {tu, tv} of parametric evaluation points. 
%                or scattered parametric coordinates (see nrleval).
%
% OUTPUT:
%
%   pnt  - interpolation matrix used to evaluate points.
%   jac  - interpolation matrix used to evaluate first derivatives (Jacobian).
%   Jacb - the determine of Jacobian at the evaluation nodes 
% 

if (nargin < 2)
  error('Not enough input arguments');
end

if (~isstruct(tsrf))
  error('NURL representation is not structure!');
end

if (~strcmp(tsrf.form,'T-NURL'))
  error('Not a recognised triangle NURL representation');
end

srf=tsrf.faces;
dd=max(srf.coefs(:))-min(srf.coefs(:));
if max(srf.coefs(3,:))-min(srf.coefs(3,:))>(1e-6)*dd
    error('Input is not a plane surface.');
end

% Coefficients (coordinates) of the triangle    
num1=srf.number(1); num2=srf.number(2); 
TN=num1*(num2-1)+1;
pnts=srf.coefs(1:3,1:TN);

% Evaluate the weights of interpolation and derivatives
[pnt, Gu, Gv]=nrltrgmat(tsrf, tt);
dpu=pnts*Gu; dpv=pnts*Gv; 

% Get first derivatives and weights of derivatives in Cartesian coordinates
Jacb=dpu(1,:).*dpv(2,:)-dpu(2,:).*dpv(1,:);
dpxi=(Gu.*repmat(dpv(2, :), [TN, 1]) ...
        - Gv.*repmat(dpu(2, :), [TN, 1]))./ ...
          repmat(Jacb, [TN, 1]); 
dpyi=(Gv.*repmat(dpu(1, :), [TN, 1]) ...
        - Gu.*repmat(dpv(1, :), [TN, 1]))./ ...
          repmat(Jacb, [TN, 1]); 

jac={dpxi, dpyi};
 

%% Demo
% % The number integration nodes (m, n)
% m=10; n=13; 
% 
% % Creat a triangular nurls patch
% crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
% af=pi/2; a=cos(af); b=sin(af);
% crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
% crv3 = nrlcirc(1, [0, 0], 0, af);
% srf = nrltrgcoons(crv1, crv2, crv3);
% 
% % Transform the surface into Lagrange blending surface
% srf = nrldegelev(srf, [2, 2]);
% srf=nrlkntins(srf, {15, 15});
% tsrf=nrl2trg(srf); 
% srf=tsrf.faces;
% 
% figure; hold on;
% nrlplot(srf, [20, 20]);
% nrlaxisplot(srf);
% view(2); axis equal;
% 
% % Get first derivatives and weights of derivatives in Cartesian coordinates
% t1=linspace(0,1,m); t2=linspace(0,1,n); 
% tt={t1, t2}; 
% [intpi, jac] = srftrgplaneder (tsrf, tt);
% dpxi=jac{1}; dpyi=jac{2}; 
% 
% % Coefficients (coordinates) of the triangle    
% num1=srf.number(1); num2=srf.number(2); 
% TN=(num2-1)*num1+1;
% pnts=srf.coefs(1:3,1:TN);
% 
% % Get a sin function and evaluate it
% x=reshape(srf.coefs(1, :), num1, num2);
% y=reshape(srf.coefs(2, :), num1, num2);
% F=sin(x).*sin(y);
% Fx=cos(x).*sin(y);
% Fy=sin(x).*cos(y);
% figure; surf(x,y,F);
% hold on; shading interp;
% Fe=sin(pnts(1,:)).*sin(pnts(2,:));
% plot3(pnts(1,:), pnts(2,:), Fe(:), 'ro'); 
% 
% xi=reshape(pnts(1,:)*intpi, m, n);
% yi=reshape(pnts(2,:)*intpi, m, n);
% Fi=reshape(Fe(:)'*intpi, m, n);
% figure; surf(xi, yi, Fi); title('Fi');
% figure; surf(xi, yi, Fi-sin(xi).*sin(yi)); 
% title('error of Fi'); 
% 
% Fxi=reshape(Fe(:)'*dpxi, m, n);
% figure; surf(xi, yi, Fxi); title('Fxi');
% figure; surf(xi, yi, Fxi-cos(xi).*sin(yi)); title('error of Fxi');
% 
% Fyi=reshape(Fe(:)'*dpyi, m, n);
% figure; surf(xi, yi, Fyi); title('Fyi');
% figure; surf(xi, yi, Fyi-sin(xi).*cos(yi)); title('error of Fyi');


