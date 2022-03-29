n = 50;
XY = 10 * rand(2,n) - 5;
for i=1:n
    plot(XY(1,i),XY(2,i),'or','MarkerSize',5,'MarkerFaceColor','r')
    axis([-5 5 -5 5])
    pause(.1)
end