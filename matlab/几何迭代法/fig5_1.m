clear 
num=40;
a=5;
theta=linspace(0,2*pi,num);
r=a*(1-cos(theta));
x=r.*cos(theta);
y=r.*sin(theta);


plot(x,y,'.--','markersize',10,'linewidth',1.1)

P=[x;y]';    

Pname=['fig5_1-',num2str(num),'.mat'];
pathname='C:\CodeStore\matlab\几何迭代法\data\';

save([pathname,Pname],'P')
axis equal