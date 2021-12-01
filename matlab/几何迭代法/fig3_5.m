clear
i=0:19;
t=i*6/19*pi;
x=3*cos(t);
y=3*sin(t).*cos(t);
z=t;

plot3(x,y,z,'.-','markersize',20);

Pname=['fig3_5-',num2str(length(i)),'.mat'];
pathname='C:\CodeStore\matlab\几何迭代法\data\';

save([pathname,Pname],'x','y','z')