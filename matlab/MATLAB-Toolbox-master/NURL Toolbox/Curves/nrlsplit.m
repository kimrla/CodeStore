function [crv1, crv2]=nrlsplit(crv, u)

% Split a nurl curve into two separate curves
% 
%  Inputs: 
% 
%      crv - the curve to be splitted
%      u  -  the knot (interval)  to split the curve
% 
%  Output: 
%      crv1 - the first half of the curve ( return the 
%                original curve is u is out of [0, 1] )
%      crv2 - the second half of the curve splitted
%                 (is empty if u is out of [0, 1] )
%
%  Examples:  
%         n=5; u=0.25;
%         crvs=nrlpolygon(n);
%         crv=nrlglues(crvs);
%         [crv1, crv2]=nrlsplit(crv, u);
%         nrlplot(crv1, 200);
%         hold on;
%         nrlplot(crv2, 200);
% 

if u<=0 || u>=1
    crv1 = crv; 
    crv2 = []; 
    return;
end

% Insert interval to preparing for partition
crv = nrlintins(crv, u);

% Get the range of each subcurve
k=find(crv.knots==u, 1, 'last');
i=find(crv.intervals==u, 1, 'last');

% Split the curve
knots1=crv.knots(1:k);
intervals1=crv.intervals(1:i);
knots1=knots1/knots1(k); 
intervals1=intervals1/intervals1(i); 
coefs1=crv.coefs(:, 1 : k);

knots2=crv.knots(k:end);
intervals2=crv.intervals(i:end);
knots2=(knots2-knots2(1))/(1-knots2(1));
intervals2=(intervals2-intervals2(1))/(1-intervals2(1)); 
coefs2=crv.coefs(:, k : end);

crv1=nrlmake(coefs1, knots1, intervals1, crv.order);
crv2=nrlmake(coefs2, knots2, intervals2, crv.order);







