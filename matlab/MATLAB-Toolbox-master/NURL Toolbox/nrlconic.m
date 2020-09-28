function [crv, d]=nrlconic(P0, T0, P2, T2, P, N)

% Create a nurl conic arc (less than pi) 
% 
%  Input:
%    P0 - start point
%    T0 - tangent vector of the point
%    P2 - end point
%    T2 - tangent vector of end point
%    P - a point between P1 and P2
%    N - the number of knots (default is 3)
% 
%  Output:
%    crv - a nurl curve
%    d - maximum distance of the lines. 
%         if d is large, the points given are 
%         not on the same plane
%

% Get points P1 and Q
u=linspace(0, 1, N);
vn=cross(T0/norm(T0), T2/norm(T2));
if norm(vn)<1e-10
    [d, Q]= planepls('neardist', P, T0, P0, P2-P0);
    p=sqrt(norm(Q-P0)/norm(P2-Q));
    u0=p/(1+p);
    f=2*u0*(1-u0)/((1-u0)^2+u0^2);
    P1=(P-Q)/f;
    w0=(1-u).^2;
    w1=2*u.*(1-u);
    w2=u.^2;    
    w=w0+w2;
else
    [d1, P1]= planepls('neardist', P0, T0, P2, T2);
    [d2, Q]= planepls('neardist', P1, P-P1, P0, P2-P0);
    d=max([d1, d2]); 
    p=sqrt(norm(Q-P0)/norm(P2-Q));
    u0=p/(1+p);
    p1=dot(P-P0, P1-P);
    p2=dot(P-P2, P1-P);
    p3=norm(P1-P)^2;
    w1=(p1*(1-u0)^2+p2*u0^2)/(2*u0*(1-u0)*p3);
    w0=(1-u).^2;
    w1=2*u.*(1-u)*w1;
    w2=u.^2;    
    w=w0+w1+w2;
end
pnts1=(P0(1)*w0+P1(1)*w1+P2(1)*w2)./w;
pnts2=(P0(2)*w0+P1(2)*w1+P2(2)*w2)./w;
pnts3=(P0(3)*w0+P1(3)*w1+P2(3)*w2)./w;
coefs=[pnts1; pnts2; pnts3; w];

crv=nrlmake(coefs, u);


%% Demo
% % Major and minor semi-axes (a, b)
% % Start and end angles (sang, eang)
% % Coordinates of the center (center)
% % The number of knots (N)
% a=2; b=1; N=12;
% sang=0; eang=2*pi/3;
% center=[0, 0];
% 
% % Preparing data for construction a nurl curve
% dang=(eang-sang);
% ang=dang*[0, 0.8, 1]+sang;
% x=a*cos(ang)+center(1);
% y=b*sin(ang)+center(2);
% dx=-dang*a*sin(ang);
% dy=dang*b*cos(ang);
% P0=[x(1), y(1), 0];
% P=[x(2), y(2), 0];
% P2=[x(3), y(3), 0];
% T0=[dx(1), dy(1), 0];
% T2=[dx(3), dy(3), 0];
% 
% % Get conic coefficients
% [crv, d]=nrlconic(P0, T0, P2, T2, P, N);
% 
% % Plot points and vectors
% ang=dang*(0:1/100:1)+sang;
% s=a*cos(ang)+center(1);
% t=b*sin(ang)+center(2);
% figure; hold on;
% plot(s, t);
% plot(x, y, 'ro')
% quiver(x, y, dx, dy);
% plot3(crv.coefs(1,:), crv.coefs(2,:), crv.coefs(3,:), 'r.');
% axis equal;
% 
% figure; hold on;
% nrlplot(crv, 1000, 'ctrl'); 
% nrlplot(crv, 5, 'quiver'); 
% axis equal;


