function [L, M, N]=quadratic2(hess, nor)

% First quadratic form of a surface
% 
% INPUTS:
%
%   hess - evaluated second derivatives (Hessian).
%   nor  -  normal vectors of the surfaces
%
% OUTPUT:
%
%   L, M, N  - Second quadratic form of a surface 
% 

L=dot(hess{1,1}(:,:), nor(:,:));
M=dot(hess{1,2}(:,:), nor(:,:));
N=dot(hess{2,2}(:,:), nor(:,:));


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
% % First quadratic form and normal vectors
% [E, F, G, g, nor]=quadratic1(jac);
% 
% % Second quadratic form
% [L, M, N]=quadratic2(hess, nor);
% 
% % Mean and Gaussian curvatures
% H=0.5*(E.*N-2*F.*M+L.*G)./g;
% K=(L.*N-M.^2)./g;





