%第三章83
clear
i=0:19;
t=3/10*pi+i/10*pi;
x=t.*cos(4*t)+3;
y=t.*sin(4*t)+3;
z=1.2*t;

plot3(x,y,z,'.-','markersize',20);

Pname=['fig3_2-',num2str(length(i)),'.mat'];
pathname='C:\CodeStore\matlab\几何迭代法\data\';

save([pathname,Pname],'x','y','z')