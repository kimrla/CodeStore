x=-2:0.05:2;
y=sin(2*x)+2*exp(-30*x.^2)+2;
dy = gradient(y);
dx = gradient(x);
% quiver(x,y,-dy,dx)
hold on; plot( x, y)
P = [x;y];
n = length(x)-1; k = 2;
NodeVector = U_quasi_uniform(n, k); % 准均匀B样条的节点矢量
DrawSpline(n, k, P, NodeVector);