function crv = nrlellip(varargin)

% Create an ellipse nurl curve
%
%   Input:
%      a - major semi-axis 
%      b - minor semi-axis
%      center - coordinates of the center
%      sang - start angle
%      eang - end angle
%
%  Output:
%      crv - an nurl curve of ellipse
%

a=varargin{1}; 
b=varargin{2};
center=varargin{3}; 
sang=varargin{4}; 
eang=varargin{5};         

% Set the angle into 2*pi
dang0=eang-sang;
dang=dang0-fix(dang0/(2*pi))*2*pi;
if abs(dang0)>0 && dang==0
    dang=2*pi;    
end
eang=sang+dang;

% Get the curves
if abs(dang)<=pi
    [pnts, u]=elliparc(a, b, sang, eang, center, 3);
    crv=nrlmake(pnts, u);
else
    mang=(sang+eang)/2;
    [pnts, u]=elliparc(a, b, sang, mang, center, 3);
    crv1=nrlmake(pnts, u);
    [pnts, u]=elliparc(a, b, mang, eang, center, 3);
    crv2=nrlmake(pnts, u);
    crv=nrlglue(crv1, crv2);
end


function [pnts, u]=elliparc(a, b, sang, eang, center, N)

% Get an elliptic arc in rational form
% 
%   Input:
%      a - major semi-axis 
%      b - minor semi-axis
%      sang - start angle
%      eang - end angle
%      center - coordinates of the center
%      N - number of knots
%   
%   Output:
%      pnts - points on the curve (first 3 rows) and weights (4th row)
%      u - knots vector defined in [0, 1]
%

% knots vector
u=linspace(0, 1, N);

% start and end points
P0=[a*cos(sang), b*sin(sang)]+center;
P2=[a*cos(eang), b*sin(eang)]+center;

% midpoint of P0-P2
M=(P0+P2)/2; 

% midpoints on the arc
mang=(sang+eang)/2;
S=[a*cos(mang), b*sin(mang)]+center;

% tangent vector at two end points
dP0=[-a*sin(sang), b*cos(sang)];
dP2=[a*sin(eang), -b*cos(eang)];

% intersection of the two tangent lines
ABC0= planeline('pointtang', P0, dP0);
ABC2= planeline('pointtang', P2, dP2);
P1=planeline('intersection', ABC0, ABC2 );

% rational curve with weights
if isnan(P1)
    w0=(1-u).^2;
    w1=2*u.*(1-u);
    w2=u.^2;
    w=w0+w2;
    P1=S-M;
else
    w1=norm(S-M)/norm(P1-S);
    w0=(1-u).^2;
    w1=2*u.*(1-u)*w1;
    w2=u.^2;
    w=w0+w1+w2;
end
pnts1=(P0(1)*w0+P1(1)*w1+P2(1)*w2)./w;
pnts2=(P0(2)*w0+P1(2)*w1+P2(2)*w2)./w;
pnts=zeros(4, N);
pnts(1,:)=pnts1(:);
pnts(2,:)=pnts2(:);
pnts(4,:)=w(:);








