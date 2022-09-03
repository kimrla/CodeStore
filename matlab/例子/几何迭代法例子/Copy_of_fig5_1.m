clear
num=500;
theta=linspace(0,2*pi,2000);
% theta=linspace(pi/2,3/2*pi,2000);
x=sin(2*theta);
y=cos(3*theta);
figure
plot(x,y,'.--','markersize',10,'linewidth',1.1)
figure
line(x,y)

P=[x;y]';
csl = [0;cumsum(vecnorm(diff(P),2,2))]; %累加长度
L=csl(end);

t=csl/L;
tt=linspace(0,1,num)';

Pp=interp1(t,P,tt,'linear');
clear P
P=Pp;
clear Pp
figure
scatter(P(:,1),P(:,2));
axis equal
% figure
% for i=1:length(P)
%     pause(0.01)
%     plot(P(i,1),P(i,2),'.','markersize',20)
%     hold on
% end
Pname=['fig5_1-',num2str(num),'.mat'];
pathname='C:\CodeStore\matlab\几何迭代法\data\';

save([pathname,Pname],'P')
