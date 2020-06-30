clear all;
close all;
clc;
N=100;%粒子个数
D=10;%粒子维度
T=200;%最大迭代次数
c1=1.5;%学习因子
c2=1.5;
w=0.8;%惯性权重
xmax=20;%位置
xmin=-20;
vmax=10;%速度
vmin=-10;
x=rand(N,D)*(xmax-xmin)+xmin;%限定位置和速度
v=rand(N,D)*(vmax-vmin)+vmin;
p=x;
pbest=ones(N,1);
for i=1:N
    pbest(i)=func1(x(i,:));
end
g= ones(1,D);
gbest = inf;
for i=1:N
    if(pbest(i)<gbest)
        g = p(i,:);
        gbest=pbest(i);
    end
end
gb=ones(1,T);
for i=1:T
    for j=1:N
        if(func1(x(j,:))<pbest(j))
            p(j,:)=x(j,:);
            pbest(j)=func1(x(j,:));
        end
        if(pbest(j)<gbest)
            g=p(j,:);
            gbest=pbest(j);
        end
        v(j,:)=w*v(j,:)+c1*rand*(p(j,:)-x(j,:))...
            +c2*rand*(g-x(j,:));
        x(j,:)=x(j,:)+v(j,:);
        for ii =1:D
            if(v(j,ii)>vmax) || (v(j,ii)<vmin)
                v(j,ii)=rand*(vmax-vmin)+vmin;
            end
            if(x(j,ii)>xmax)||(x(j,ii)<xmin)
                x(j,ii)=rand*(xmax-xmin)+xmin;
            end
        end
    end
    gb(i)=gbest;
end
    g;
    gb(end);
    figure
    plot(gb)
    xlabel('迭代次数');
    ylabel('适应度');
    title('适应度进化曲线')