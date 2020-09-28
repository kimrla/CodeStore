function [Ke, Me, pnts, G, Gx, Gy]=StiffMassInpTri(tsrf, E, nu, rho)

% StiffnessMass: Get stiffness and mass matrices for isotropic plane problem on triangles
% 
% Calling Sequences:
%
%     [Ke, Me, pnt, G, Gx, Gy]=StiffMassInpTri(srf, E, nu, rho)
%
% INPUTS:
%
%      tsrf - a nurl surface for triangle
% 
%      E, nu, rho -  Young's modulus (E), Poisson's ratio (nu) and density (rho)
%
% OUTPUT:
%
%   Ke, Me  -  stiffness and mass matrices
%   pnts  - points on integration nodes
%   G, Gx, Gy - weighting coefficient matrix for evaluation and derivatives
% 

% Get first derivatives and weights of derivatives in Cartesian coordinates
[G, Gx, Gy, ~, ~, C, pnts, Jacb]=nrlplanemattrg(tsrf);

% Stiffness and mass matrices
CJ=diag(C(:).*Jacb(:));
Cr=E/(1-nu^2); nu1=(1-nu)/2;
G11=Gx'*CJ*Gx;
G22=Gy'*CJ*Gy;
G12=Gx'*CJ*Gy;
Ke=Cr*[G11+nu1*G22, nu*G12+nu1*G12';
           nu*G12'+nu1*G12, G22+nu1*G11];
Gm=G'*CJ*G; 
Zm=zeros(size(Gm)); 
Me=rho*[Gm, Zm; Zm, Gm]; 







