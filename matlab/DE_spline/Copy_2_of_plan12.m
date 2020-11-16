clear all;

NP=50;%种群规模NP

GM=200;%最大迭代次数

p=3;%B样条次数p=3，控制顶点n+1个，节点矢量ui i=1~n+p+2,

% % 实验1
% x=0:0.005:1-0.005;
% f=90./(1+exp(-100*(x-0.4)));
% n=8;%内节点个数n-p
% lamda=0.025;%节点率λ
% dp=0.6;%删除概率
% Num=length(x);%采样点个数
% a=min(x);
% b=max(x);

% 实验2
x=0:0.05:10-0.05;
f=100./exp(abs(x-5))+(x-5).^2/500;
n=13;
lamda=0.05;
dp=0.6;
Num=length(x);
a=min(x);
b=max(x);

f_=f+normrnd(0,0.01,1,Num);%加随机扰动


% plot(x,f_,"*")
% hold on
d=[x;f_]';%给定数据点,di=(xi,yi) i=1~M+1


M=length(x)-1; %论文中是0~M，所以总数length=M+1


p=3;%B样条次数p=3，控制顶点n+1个，节点矢量ui i=1~n+p+2,

%初始种群


for i=1:NP
    X{i}=sort(neijiedian(n,p,a,b)+normrnd(0,0.01,1,n-3));%初始种群并初始化
    
    ui{i}=[zeros(1,p+1) X{i} b*ones(1,p+1)];%方案1 均匀节点向量     
    [N{i},~,P{i}] = kongzhidingdian(M,n,p,x,ui{i},d);
    [~,R(i)]=shujudianwucha(M,N{i},P{i},d);
    BIC(i)=Num*log(1+R(i))+25*log(Num)*(2*n-p+1);     
end
ui_=cell (1,50);
for gen=1:GM
    for m=1:NP       
        r1=randi([1,NP],1,1);
        while(r1==m)
            r1=randi([1,NP],1,1);
        end
        
        r2=randi([1,NP],1,1);
        while(r2==r1)||(r2==m)
            r2=randi([1,NP],1,1);
        end
        
        r3=randi([1,NP],1,1);
        while(r3==m)||(r3==r2)||(r3==r1)
            r3=randi([1,NP],1,1);
        end
        r=[r1,r2,r3];
        %  产生不同的r1,r2,r3
        for i=1:3
            rand1=rand;%增减节点数
            X_{i}=X{r(i)};
            if rand1<dp
                X_{i}(randi(length(X_{i}),1,1))=[];
%             else                
%                 X_{i}(end+1)=a+(b-a)*rand;
%                 X_{i}=sort(X_{i});
            end
            nr(i)=length(X_{i});
        end
       D=min(nr);
       for i=1:3
           for k=1:D
               temp{i}(k)=sum(X_{i}(k:k+nr(i)-D))/(nr(i)-D+1);    
           end
           X_{i}=temp{i}; 
       end
       temp=cell(1,3);
       
       FG=0.5*exp(GM/(GM+gen)-1);
       U{m}=X_{1}+(X_{2}-X_{3})*FG;
       
       %杂交操作
       
       for j=1:D
           rand4=rand;
           rand2=rand;
           rand3=rand;
           CR=0.5*(rand4+1);
           if j<=length(X{m}) & rand2>CR 
               V{m}(j)=X{m}(j);
           elseif  j<=D & U{m}(j)>a & U{m}(j)<b
               V{m}(j)=U{m}(j);
           else
               V{m}(j)=a+(b-a)*rand3;
           end
       end
       V{m}=sort(V{m});
       
       
       ui_{m}=[zeros(1,p+1) V{m} b*ones(1,p+1)];%方案1 均匀节点向量       
       n_=length(ui_{m})-p-2;
       [N_{m},~,P_{m}] = kongzhidingdian(M,n_,p,x,ui_{m},d);
       [~,R_(m)]=shujudianwucha(M,N_{m},P_{m},d);
       BIC_(m)=Num*log(1+R_(m))+25*log(Num)*(2*n_-p+1);     
       
       if BIC_(m)<BIC(m)
           X{m}=V{m};
           BIC(m)=BIC_(m);
           P{m}=P_{m};
           R(m)=R_(m);
           ui{m}=ui_{m};
       end
    end
    [BICbest(gen),position(gen)]=min(BIC);
    neijiedianshuliang(gen)=length(X{position(gen)});
    bestP{gen}=P{position(gen)};
    bestui{gen}=ui{position(gen)};
    
end
figure
plot(x,f_,"*")
hold on
DrawSpline(neijiedianshuliang(end)+p,p,bestP{end},bestui{end},a,b);
scatter(bestui{end},ones(1,length(bestui{end})))
figure
plot(BICbest)
xlabel('迭代次数');
ylabel('适应度');
title('适应度进化曲线')
figure
plot(neijiedianshuliang)



