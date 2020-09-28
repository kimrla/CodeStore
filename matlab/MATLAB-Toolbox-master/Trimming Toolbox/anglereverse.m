function theta=anglereverse(pt21)
% Reverse the theta. The range of the angle can be [0,2*pi), compared with
% the MATLAB function asin (whose range is [-pi/2,pi/2]) or acos (whose
% range is [0,pi]).

% Input:
%     pt21: The direction vector, which is a 1*2 vector for 2D.
% Output:
%     theta: Angle from x-axis in global co sys to pt21.

pt21=vecnorm(pt21);
costheta=dot([1,0,0]',[pt21(:)',0]);
sintheta=cross([1,0,0]',[pt21(:)',0]);
sintheta=sintheta(3);
theta1=acos(costheta);
theta2=asin(sintheta);
if abs(theta1-theta2)<=10e-6 % The first quadrant
    theta=theta1;
else
    if theta1>=pi/2 && theta2>0 % The second quadrant
        theta=theta1;
    elseif theta1>=pi/2 && theta2<=0 % The third quadrant
        theta=pi-theta2;
    elseif theta1<=pi/2 && theta2<0 % The fourth quadrant
        theta=2*pi+theta2;
    end
end
    
end