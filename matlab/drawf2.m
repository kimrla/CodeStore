x=0:0.4:16*pi;
y=sin(x);
dy = gradient(y);
dx = gradient(x);
% quiver(x,y,-dy,dx,'autoscale','on')
hold on; plot( x, y)