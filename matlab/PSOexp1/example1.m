clear all;
close all;
clc;
N=100;%���Ӹ���
D=10;%����ά��
T=200;%����������
c1=1.5;%ѧϰ����
c2=1.5;
w=0.8;%����Ȩ��
xmax=20;%λ��
xmin=-20;
vmax=10;%�ٶ�
vmin=-10;
x=rand(N,D)*(xmax-xmin)+xmin;%�޶�λ�ú��ٶ�
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
    xlabel('��������');
    ylabel('��Ӧ��');
    title('��Ӧ�Ƚ�������')