x=-2:0.04:2;
y=sin(2*x)+2*exp(-30*x.^2)+2;
dy = gradient(y);
dx = gradient(x);
quiver(x,y,-dy,dx)
hold on; plot( x, y)