%几何迭代法第二章37页
clear
a=1;
num=10;%采样点个数
theta=linspace(0,2*pi,num);
x=a*(1+cos(theta));
y=a*sin(theta).*(1+cos(theta));
plot(x,y,'.-','markersize',20)

P=[x;y]';    

Pname=['fig2_3-',num2str(num),'.mat'];
pathname='C:\CodeStore\matlab\几何迭代法\data\';

save([pathname,Pname],'P')
