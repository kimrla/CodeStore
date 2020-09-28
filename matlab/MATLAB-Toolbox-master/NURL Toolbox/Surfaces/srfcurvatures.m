function [H, K]=srfcurvatures(jac, hess)

% Mean and Gaussian curvatures of a surface
% 
% INPUTS:
%
%   jac  - evaluated first derivatives (Jacobian).
%   hess - evaluated second derivatives (Hessian).
% 
% OUTPUT:
%
%   H  - Mean curvatures
%   K  - Gaussian curvatures
% 

% First quadratic form and normal vectors
[E, F, G, g, nor]=quadratic1(jac);

% Second quadratic form
[L, M, N]=quadratic2(hess, nor);

% Mean and Gaussian curvatures
H=0.5*(E.*N-2*F.*M+L.*G)./g;
K=(L.*N-M.^2)./g;


%% Demo
% R=1; 
% s1=0; s2=pi/2; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;
% 
% % First and second derivatives
% m=10; n=11;
% t1=linspace(0, 1, m);
% t2=linspace(0, 1, n);
% tt={t1, t2};
% [pnts, jac, hess]=nrldeval(srf, tt);
% 
% % Mean and Gaussian curvatures
% [H, K]=srfcurvatures(jac, hess);
% 
% % Principal curvatures 
% dt=real(sqrt(H.^2-K));
% k1=H+dt;
% k2=H-dt;


