function [a, ax, ay]=plangle(v)

% Get the angle of a plane vector with respect to x-coordinate
%    input       :  v - vector
%    output    :  a - angle
%    example :  v=[ 0.5, 1]; 
%                     a =angle(v);

vx=[1, 0]; vy=[0, 1];

% Function handles used to compute angles
radvec = @(v) sqrt(sum(sum(v.^2))); 
vtov = @(v1, v2) acos(dot(v1,v2)/radvec(v1)/radvec(v2))*180/pi; 

% Get the angle with respecto vx and vy
ax = vtov(v, vx); 
ay = vtov(v, vy); 

% Judge the quadrant of the vector
if ay>=0 && ay<90
    a=ax;
else
    a=360-ax;
end

a=a*pi/180;

