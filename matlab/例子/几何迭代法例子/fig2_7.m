%第二章46页
clear
num=20;%采样点个数
theta=linspace(-3*pi,3*pi,num);
x=cos(theta);
y=sin(theta);
z=theta/2;
plot3(x,y,z,'.-','markersize',20);

Pname=['fig2_7-',num2str(num),'.mat'];
pathname='C:\CodeStore\matlab\几何迭代法\data\';

save([pathname,Pname],'x','y','z')