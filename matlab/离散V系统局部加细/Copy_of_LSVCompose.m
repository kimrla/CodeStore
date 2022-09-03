function LSVCompose(Lambda,k,N)
%% V基合成 Lambda：系数；N：组数
 % 方法：在每段子区间上采样k+1个点，用k次Lagrange多项式插值得到该子区间上多项式

%%
[NumLda,NumDim] = size(Lambda);  % 给定系数长度及维度
NumV = (k+1)*2^(N-1);            % 合成基函数数目

% if NumLda ~= NumV
%     error('ERROR!!!The Number of coefficient is not match with basis!');
% end

NumSeg = 2^(N-1);                % 分段数
NumSmp = (k+1)*NumSeg;           % 每个子区间采样k+1个点
DltSmp = 1/NumSmp;               % 采样间隔
T = [DltSmp/2 : DltSmp : 1]';    % 采样节点

%% 合成
A = LSMatrix_V(k,N,T); 
Vdiscrete = A(:,1:NumLda) * Lambda;     

%% 逐段插值为k次多项式
NumSmpOnSeg = 65;                                 % 每段采样点数（默认）
LengthSubinterval = 1/NumSeg;
t0 = linspace(0,LengthSubinterval,NumSmpOnSeg);   % 第1段采样节点

for s = 1 : NumSeg
    ts = t0 + (s-1)*LengthSubinterval;
    ys = zeros(NumSmpOnSeg,NumDim);
    Vs = Vdiscrete((s-1)*(k+1)+1:s*(k+1),:);
    for dim = 1 : NumDim
        ys(:,dim) = lagrangepoly(T((s-1)*(k+1)+1:s*(k+1)),Vs(:,dim),ts);
    end 
    switch NumDim
        case 1
            plot(ts,ys,'Color',[0 102 153]/255,'LineWidth',1.5);hold on
        case 2
            plot(ys(:,1),ys(:,2),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
        case 3
            plot(ys(:,1),ys(:,2),ys(:,3),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
    end
end


function [P,R,S] = lagrangepoly(X,Y,XX)
%LAGRANGEPOLY  Lagrange interpolation polynomial fitting a set of points
%   [P,R,S] = LAGRANGEPOLY(X,Y)  where X and Y are row vectors
%   defining a set of N points uses Lagrange's method to find 
%   the N-1th order polynomial in X that passes through these 
%   points.  P returns the N coefficients defining the polynomial, 
%   in the same order as used by POLY and POLYVAL (highest order first).
%   Then, polyval(P,X) = Y.  R returns the x-coordinates of the N-1
%   extrema of the resulting polynomial (roots of its derivative),
%   and S returns the y-values  at those extrema.
%
%   YY = LAGRANGEPOLY(X,Y,XX) returns the values of the polynomial
%   sampled at the points specified in XX -- the same as
%   YY = POLYVAL(LAGRANGEPOLY(X,Y)). 
%
%   Example:
%   To find the 4th-degree polynomial that oscillates between 
%   1 and 0 across 5 points around zero, then plot the interpolation
%   on a denser grid inbetween:
%     X = -2:2;  Y = [1 0 1 0 1];
%     P = lagrangepoly(X,Y);
%     xx = -2.5:.01:2.5;
%     plot(xx,polyval(P,xx),X,Y,'or');
%     grid;
%   Or simply:
%     plot(xx,lagrangepoly(X,Y,xx));
%
%   Note: if you are just looking for a smooth curve passing through 
%   a set of points, you can get a better fit with SPLINE, which 
%   fits piecewise polynomials rather than a single polynomial.
%
%   See also: POLY, POLYVAL, SPLINE

% 2006-11-20 Dan Ellis dpwe@ee.columbia.edu
% $Header: $

%  For more info on Lagrange interpolation, see Mathworld: 
%  http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html

% Make sure that X and Y are row vectors
if size(X,1) > 1;  X = X'; end
if size(Y,1) > 1;  Y = Y'; end
if size(X,1) > 1 || size(Y,1) > 1 || size(X,2) ~= size(Y,2)
  error('both inputs must be equal-length vectors')
end

NN = length(X);

pvals = zeros(NN,NN);

% Calculate the polynomial weights for each order
for i = 1:NN
  % the polynomial whose roots are all the values of X except this one
  pp = poly(X( (1:NN) ~= i));
  % scale so its value is exactly 1 at this X point (and zero
  % at others, of course)
  pvals(i,:) = pp ./ polyval(pp, X(i));
end

% Each row gives the polynomial that is 1 at the corresponding X 
% point and zero everywhere else, so weighting each row by the 
% desired row and summing (in this case the polycoeffs) gives 
% the final polynomial
P = Y*pvals;

if nargin==3
  % output is YY corresponding to input XX
  YY = polyval(P,XX);
  % assign to output
  P = YY;
end

if nargout > 1
  % Extra return arguments are values where dy/dx is zero
  % Solve for x s.t. dy/dx is zero i.e. roots of derivative polynomial
  % derivative of polynomial P scales each power by its power, downshifts
  R = roots( ((NN-1):-1:1) .* P(1:(NN-1)) );
  if nargout > 2
    % calculate the actual values at the points of zero derivative
    S = polyval(P,R);
  end
end
end




end