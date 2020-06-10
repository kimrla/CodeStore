x = 0:0.1:1;
y = x.*x;
dy = gradient(y);
dx = gradient(x);
quiver(x,y,-dy,dx)
hold on; plot( x, y)