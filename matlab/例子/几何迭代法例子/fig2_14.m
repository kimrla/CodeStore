%第二章73页
clear
i=0:10;
t=-pi/2+i*pi/5;
x=cos(t);
y=sin(t).*cos(t);
plot(x,y,'.-','markersize',20)
P=[x;y]';    

Pname=['fig2_14-',num2str(length(i)),'.mat'];
pathname='C:\CodeStore\matlab\几何迭代法\data\';

save([pathname,Pname],'P')
