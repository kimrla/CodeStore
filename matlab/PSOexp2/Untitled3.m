clear all;
close all;
clc;
n=100;
d=2;
t=200;
c1=1.5;
c2=1.5;
wmax=0.8;
wmin=0.4;
xmax=4;
xmin=-4;
vmax=1;
vmin=-1;
x=rand(n,d)*(xmax-xmin)+xmin;
v=rand(n,d)*(vmax-vmin)+vmin;
p=x;
pbest=ones(n,1);
a=[];
b=[];
c=[];
for i=1:n
    pbest(i)=func2(x(i,:));
end
g=ones(1,d);
gbest=inf;
for i=1:n
    if(pbest(i)<gbest)
        g=p(i,:);
        gbest=pbest(i);
    end
end
gb=ones(1,t);
for i=1:t    
    for j=1:n
        a(j,i)=p(j,1);
        b(j,i)=p(j,2);
        c(j,i)=pbest(j);
        if(func2(x(j,:))<pbest(j))
            p(j,:)=x(j,:);
            pbest(j)=func2(x(j,:));
        end
        if(pbest(j)<gbest)
            g=p(j,:);
            gbest=pbest(j);
        end
        w=wmax-(wmax-wmin)*i/t;
        v(j,:)=w*v(j,:)+c1*rand*(p(j,:)-x(j,:))...
            +c2*rand*(g-x(j,:));
        x(j,:)=x(j,:)+v(j,:);
        for ii=1:d
            if(v(j,ii)>vmax)||(v(j,ii)<vmin)
                v(j,ii)=rand*(vmax-vmin)+vmin;
            end
            if(x(j,ii)>xmax)||(x(j,ii)<xmin)
                x(j,ii)=rand*(xmax-xmin)+xmin;
            end
        end       
    end
    gb(i)=gbest;
%     plot3(a(:,i),b(:,i),c(:,i),'.'),pause(0.01);
end
% h=plot3(a,b,c,'r');
% for m = 2:numel(a)
%    set(h,'XData',a(1:m),'YData',b(1:m),'ZData',c(1:m));
%    drawnow;
% end
% x0=-4:0.08:4;
% y0=-4:0.08:4;
% for i=1:n
%     for j=1:n
%         z0(i,j)=3*cos(x0(i)*y0(j))+x0(i)+y0(j)^2;
%         
%     end
% end
% mesh(x0,y0,z0);
% f=@(x0,y0) 3*cos(x0*y0)+x0+y0^2;
% fsurf(f,[-4 4 -4 4]);
% hold on;
% [x0,y0]=meshgrid(-4:.1:4);
% z0=3*cos(x0*y0)+x0+y0.^2;
% mesh(x0,y0,z0);
% hold on;
for i=1:t
      
%     z0=-10:0.3:20;
%     set(gca,'xlim',[-5 5]);
%     set(gca,'ylim',[-5 5]);
%     set(gca,'zlim',[-20 20]);
%     plot3(x0,y0,z0,'w',a(:,i),b(:,i),c(:,i),'.'),pause(0.01);
    plot3(a(:,i),b(:,i),c(:,i),'.'),pause(0.01);
%     for m = 2:numel(a)
%         set(h,'XData',a(1:m),'YData',b(1:m),'ZData',c(1:m));
%         drawnow;
%     end
    hold on;
end

g;
gb(end);
figure
plot(gb)
xlabel('迭代次数');
ylabel('适应度');
title('适应度进化曲线');