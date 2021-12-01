clear 
num=100001;
theta=linspace(0,40*pi,num);
r=theta;
x=r.*cos(theta);
y=r.*sin(theta);


plot(x,y,'linewidth',1.1)

P=[x;y]';    

Pname=['fig4_8-',num2str(num),'.mat'];
pathname='C:\CodeStore\matlab\几何迭代法\data\';

save([pathname,Pname],'P')
