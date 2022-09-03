clear 
num=501;
theta=linspace(0,8*pi,num);
r=sin(theta/4);
x=r.*cos(theta);
y=r.*sin(theta);


plot(x,y,'.-','markersize',20)

P=[x;y]';    

Pname=['fig4_7-',num2str(num),'.mat'];
pathname='C:\CodeStore\matlab\几何迭代法\data\';

save([pathname,Pname],'P')
