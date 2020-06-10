t=-4*pi:0.2:4*pi;
x=sin(0.75*t);
y=sin(t);
dy = gradient(y);
dx = gradient(x);
quiver(x,y,-dy,dx)
hold on; plot( x, y)