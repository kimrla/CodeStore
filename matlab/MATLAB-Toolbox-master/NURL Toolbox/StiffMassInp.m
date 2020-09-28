function [Ke, Me, pnt, G, Gx, Gy]=StiffMassInp(srf, E, nu, rho)

% StiffnessMass: Get stiffness and mass matrices for isotropic plane problem on rectangles
% 
% Calling Sequences:
%
%     [Ke, Me, pnt, G, Gx, Gy]=StiffMassInp(srf, E, nu, rho)
%
% INPUTS:
%
%      srf - a nurl surface
% 
%      E, nu, rho -  Young's modulus (E), Poisson's ratio (nu) and density (rho)
%
% OUTPUT:
%
%   Ke, Me  -  stiffness and mass matrices
%   pnt  - points on integration nodes
%   G, Gx, Gy - weighting coefficient matrix for evaluation and derivatives
% 

% Get first derivatives and weights of derivatives in Cartesian coordinates
[s, Cs]=GaussLobattoQ([0,1], srf.order(1)+2);
[t, Ct]=GaussLobattoQ([0,1], srf.order(2)+2);
[Ct, Cs]=meshgrid(Ct, Cs);
[intpi, jac, geom] = nrlsrfgintvdeval (srf.order, srf.coefs, srf.knots, {s, t});
G=intpi'; Gx=jac{1}'; Gy=jac{2}'; 
CJ=geom{1}.*Cs(:)'.*Ct(:)';
CJ=repmat(CJ', 1, size(intpi, 1));

% Stiffness and mass matrices
C=E/(1-nu^2); nu1=(1-nu)/2;
G11=Gx'*(CJ.*Gx);
G22=Gy'*(CJ.*Gy);
G12=Gx'*(CJ.*Gy);
Ke=C*[G11+nu1*G22, nu*G12+nu1*G12';
           nu*G12'+nu1*G12, G22+nu1*G11];
Gm=G'*(CJ.*G); 
Zm=zeros(size(Gm)); 
Me=rho*[Gm, Zm; Zm, Gm]; 
pnt=nrleval(srf, {s, t});



