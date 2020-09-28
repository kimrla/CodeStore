function [crvs, u0]=nrlsplits(crv, u)

% Split a nurl curve into several separate curves
% 
%  Inputs: 
% 
%      crv - the curve to be splitted
% 
%      u - the knots to split the curve
% 要分割的节点在(0,1)之间且没有重复
%  Output: 
% 
%      crvs - the curves splitted from crv
% 
%      uu - the used konts to split the curve
% uu是筛选后符合要求的u
%  Examples:  
% 
%         n=5; u=[-0.1, 0.375,0.25,0.575,0.6375,0.85, 1.1];
%         crvs=nrlpolygon(n);
%         crv=nrlglues(crvs);
%         nrlplot(crv,1000);
%         [crvs, uu]=nrlsplits(crv, u);
%         figure; hold on;
%         for i=1:numel(crvs)
%             nrlplot(crvs(i), 200);
%         end
% 

% Remove knots out of span
qq=(u<=0) | (u>=1);
u(qq)=[];

% Sort and remove repeated knots
u=sort(u(:)); 
u=RemDuplicate(u);
u=u'; u0=u;

% Split the curve into (n+1) curves        
n=length(u); 
if n>0
    [crv1, crv2]=nrlsplit(crv,u(1));
    crvs=[crv1, crv2];
    u=(u-u(1))/(1-u(1));
    for i=2:n
        [crv1, crv2]=nrlsplit(crvs(i),u(i));
        crvs(i)=[];
        crvs=[crvs, crv1, crv2];
        u=(u-u(i))/(1-u(i));
    end
else
    crvs=crv;
end







