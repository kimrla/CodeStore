function crv=nrlpolygon(n, radius, center, ang)

% Create a nurl polygon in x-y plane
% 
%   Inputs : 
% 
%      radius, center, ang
%      radius : Radius of the circle, default 1.0
%      center : Center of the circle, default (0,0)
%      ang : Angle of the first line with respect x-coordinate, default 0 radians (0 degrees)
% 
%  Example   :  
% 
%     n=5; radius=1; center=[0,0]; ang=0;
%     crv=nrlpolygon(n,radius,center,ang);
%     crv=nrlpolygon(n,radius,center);
%     crv=nrlpolygon(n,radius);
%     crv=nrlpolygon(n);
% 

m=nargin;
if m==1
    radius=1.0;
    center=[0,0];
    ang=0;
elseif m==2
    center=[0,0];
    ang=0;
elseif m==3
    ang=0;
end

deta=pi/n;
xy=zeros(n,2);
for i=1:n
    angi=ang+2*i*deta-deta;
    xy(i,1)=center(1)+radius*sin(angi);
    xy(i,2)=center(2)-radius*cos(angi);
end
crv(1,n)=nrlline(xy(n,:), xy(1,:));
for i=1:n-1
    crv(i)=nrlline(xy(i,:), xy(i+1,:));
end    



%% Demo
% n=5; radius=1; center=[0,0]; ang=0;
% crv=nrlpolygon(n, radius, center, ang);
% figure; hold on;
% for i=1:n
%     nrlplot(crv(i));
% end





