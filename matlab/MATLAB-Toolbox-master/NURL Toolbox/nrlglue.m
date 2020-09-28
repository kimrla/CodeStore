function [crv, id]=nrlglue(crv1, crv2)

% Glue two nurl curves together
% 
%  Input:
%    crv1 - curve 1
%    crv2 - curve 2
%
%  Output:
%    crv - the curve glued crv1 and crv2 together
%

% Get index of the common points of two curves
stedpnts1=nrleval(crv1, [0, 1]);
stedpnts2=nrleval(crv2, [0, 1]);
DM = DistanceMatrix(stedpnts1', stedpnts2');
dist = nrlmeasure (crv1);
index=abs(DM)<dist*1e-3;
ii=[1,1;2,2]; jj=[1,2;1,2];
i=ii(index); j=jj(index);
if numel(i)==0
    crv = []; 
    id = 0; 
    return;
end
i=i(1); j=j(1);

% Make the end of first curve be the start of the second curve
if i==1 && j==1
    crv1 = nrlreverse(crv1);
elseif i==1 && j==2
    crv1 = nrlreverse(crv1);
    crv2 = nrlreverse(crv2);
elseif i==2 && j==2
    crv2 = nrlreverse(crv2);
end

% Elevate the degrees
p=crv1.order-crv2.order;
if p>0
    crv2=nrldegelev(crv2, p);
elseif p<0
    crv1=nrldegelev(crv1, -p);
end

% Glue the curves
len1=nrlmeasure(crv1);
len2=nrlmeasure(crv2);
d=len1/(len1+len2);
n=length(crv2.knots);
m=length(crv2.intervals);
knots=[d*crv1.knots, (1-d)*crv2.knots(2:n)+d];
coefs=[crv1.coefs, crv2.coefs(:, 2:n)];
intervals=[d*crv1.intervals, (1-d)*crv2.intervals(2:m)+d];
k=length(crv1.intervals);
p=length(crv1.knots);
knots(p)=intervals(k);
crv=nrlmake(coefs, knots, intervals, crv1.order);
id=1;



