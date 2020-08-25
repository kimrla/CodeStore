clear all;
% 实验1
x=-2:0.05:2-0.05;
y=sin(2*x)+2*exp(-30*x.^2)+2;
% % 实验2
% x=0:0.4:16*pi;
% y=sin(x);
% 实验3
% t=-4*pi:0.2:4*pi;
% x=sin(0.75*t);
% y=sin(t);


d=[x;y]';%给定数据点,di=(xi,yi) i=1~M+1

dy = gradient(y);
dx = gradient(x);
l=[-dy;dx]';%在数据点处的法向量约束条件li i=1~M+1

M=length(x)-1; %论文中是0~M，所以总数length=M+1

t=canshuhua(M,d);
p=3;%B样条次数p=3，控制顶点n+1个，节点矢量ui i=1~n+p+2  n=M+1+p-1,

n=30;

P=60;%粒子规模P 内节点数量
T=10;%最大迭代次数
c1=2;
c2=2;
w=0.7;
for i=1:P
    X(i,:)=sort(rand(1,n-p));%初始位置，内节点
end
% V=sort(rand(P,n-p));%初始速度
V=zeros(P,n-p);
xmax=1;%位置
xmin=0;
ui=[zeros(P,p+1) X ones(P,p+1)];%PSO自适应节点向量
% ui=jiedianxiangliang(n,p);%方案1 均匀节点向量
% ui=jiedianxiangliang2(n,p,t);%方案2 平均节点向量
% ui=jiedianxiangliang3(n,p,t);%方案3 皮格尔逼近节点向量


% 最小二乘法求控制顶点
for k=1:P
    [N(:,:,k),R(:,:,k),Pi(:,:,k)] = kongzhidingdian(M,n,p,t,ui(k,:),d);


% DrawSpline(n,p,Pi,ui);
% axis([x(1)-1,x(end)+1,min(y)-1,max(y)+1])

[epsilon(:,:,k),e(k)] = shujudianwucha(M,N(:,:,k),Pi(:,:,k),d);
[sst(:,:,k),es(k,:),ev(k)] = faxiangliangwucha(M,n,p,ui(k,:),t,l,Pi(:,:,k));
end
position=X;
pbest=ones(P,1);
for i=1:P
    pbest(i)=func(M,n,p,ui(i,:),t,l,d);
end
g= ones(1,n-p);
gbest = inf;
% for i=1:P
%     if(pbest(i)<gbest)
%         g = position(i,:);
%         gbest=pbest(i);
%     end
% end
[gbest,i]=min(pbest);
g = position(i,:);
gb=ones(1,T);
for i=1:T
            if i>1
               [B,index]= sort(fit);
               X(index(P/2+1:end),:)=X(index(1:P/2),:);
            end
    for j=1:P
       
        ui(j,:)=[zeros(1,p+1) X(j,:) ones(1,p+1)];
        fit(j)=func(M,n,p,ui(j,:),t,l,d);
        
        if fit(j)<pbest(j)
            position(j,:)=X(j,:);
            pbest(j)=fit(j);
        end
        if(pbest(j)<gbest)
            g=position(j,:);
            gbest=pbest(j);
        end
        V(j,:)=w*V(j,:)+c1*rand*(position(j,:)-X(j,:))...
            +c2*rand*(g-X(j,:));
        X(j,:)=X(j,:)+V(j,:);
        for ii =1:n-p
%             if(V(j,ii)>vmax) || (V(j,ii)<vmin)
%                 V(j,ii)=rand*(vmax-vmin)+vmin;
%             end
            if(X(j,ii)>xmax) || (X(j,ii)<xmin)
                
                if ii==1
                    for iii=2:n-p
                        if (X(j,iii)<xmax) && (X(j,iii)>xmin)
                            X(j,1:iii-1)=sort(rand(1,iii-1))*(X(j,iii)-0)+0;
                        end
                    end                        
                elseif ii==n-p
                    X(j,ii)=rand*(1-X(j,ii-1))+X(j,ii-1);
                else
                    for iii=ii+1:n-p
                        if (X(j,iii)<xmax) && (X(j,iii)>xmin)
                            X(j,ii:iii-1)=sort(rand(1,iii-ii))*(X(j,iii)-X(j,ii-1))+X(j,ii-1);
                        end
                    end                
                end
                
            end
            
        end
    end
    gb(i)=gbest;
end
    g;%最优个体
    gb(end);%最优值
    figure
    plot(gb)
    xlabel('迭代次数');
    ylabel('适应度');
    title('适应度进化曲线')
fprintf('数据点处误差是%d\n法向量误差是%d',e,ev)