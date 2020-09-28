function [G, Gx, Gy, S, T, C, pnt, Jacb]=nrlplanemattrg(tsrf)

% Get first derivatives and weights of derivatives in Cartesian coordinates
% 
%  INPUT:
%
%    tsrf    -   a nurl triangle patch.
%
% OUTPUT:
% 
%    G - a matrix of the basis
% 
%    Gx, Gy - matrices of first derivatives with respect to Cartesian coordinates
% 
%    S, T - area coordinates transformed from natural coodinates tt
%
%    C - weights of integration
%
%    pnt - points on the triangular surface
%
%    Jacb -  the determinant |J| of the Jacobian J = d(x, y)/d(s, t)
%

% The total number of nodes
srf=tsrf.faces;
num1=srf.number(1); 
num2=srf.number(2); 
TN=num1*(num2-1)+1;

% Get first derivatives and weights of derivatives in Cartesian coordinates
[G, Gu, Gv, S, T, C, s, t]=nrlmattrg(tsrf); 

% Get first derivatives and weights of derivatives in Cartesian coordinates
[pnt, jac]=nrltrgdeval(tsrf, {s, t});
dpu=jac{1}; dpv=jac{2};
Jacb=dpu(1,:).*dpv(2,:)-dpu(2,:).*dpv(1,:);
dpxi=(Gu'.*repmat(dpv(2, :), [TN, 1]) ...
        - Gv'.*repmat(dpu(2, :), [TN, 1]))./ ...
          repmat(Jacb, [TN, 1]); 
dpyi=(Gv'.*repmat(dpu(1, :), [TN, 1]) ...
        - Gu'.*repmat(dpv(1, :), [TN, 1]))./ ...
          repmat(Jacb, [TN, 1]); 

Gx=dpxi';
Gy=dpyi';






