function [dr, q, dr1, dr2]=BisectVector(dr1, dr2)

% Get the bisector vector of two plane vectors
% 
%  Input : 
%
%    dr1 - vector1
% 
%    dr2 - vector 2
% 
%  Output : 
%
%    dr - bisector vector
%    q - bisector angle of the two vectors
%    dr1 - normalized dr1
%    dr2 - normalized dr2
% 
%  Example: 
%  
%     dr1=[ 0; 1]; dr2=[-1; 0];
%     [dr, q, dr1, dr2]=BisectVector(dr1, dr2);
% 

dr1=vecnorm(dr1');
dr2=vecnorm(dr2');
if norm(dr1+dr2)==0
    dr(1)=-dr1(2);
    dr(2)=dr1(1);
    q=pi/2;
else
    dr=vecnorm(dr1+dr2);
     q=vtov(dr1, dr);
end

% Angle between vector v1 and v2
function y=vtov(v1,v2,flag)
if nargin==2
    flag=1;
end
if flag==1
    y=acos(dot(v1,v2)/radvec(v1)/radvec(v2));
elseif flag==2
    y=acos(dot(v1,v2)/radvec(v1)/radvec(v2))*180/pi;
end

% Length of the vector
function y=radvec(v)
y=sqrt(sum(sum(v.^2)));




