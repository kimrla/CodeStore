clear all;

NP=50;%种群规模NP

GM=200;%最大迭代次数

% 实验1
x=0:0.01:1-0.01;
f=90./(1+exp(-100*(x-0.4)));
n=8;
lamda=0.025;%节点率λ
dp=0.6;%删除概率
N=length(x);%采样点个数

% % 实验2
% x=0:0.05:10;
% y=100./exp(abs(x-5))+(x-5).^2/500;
% n=13;

f_=f+normrnd(0,0.01,1,N);%加随机扰动


plot(x,f_,"*")
hold on
d=[x;f_]';%给定数据点,di=(xi,yi) i=1~M+1


M=length(x)-1; %论文中是0~M，所以总数length=M+1

t=canshuhua(M,d);
p=3;%B样条次数p=3，控制顶点n+1个，节点矢量ui i=1~n+p+2,

%初始种群
ui=zeros(NP,n+p+2);
for i=1:NP
    ui(i,:)=sort(jiedianxiangliang(n,p)+[zeros(1,4) normrnd(0,0.01,1,n-3) zeros(1,4)]);%方案1 均匀节点向量
end


[N,R,P] = kongzhidingdian(M,n,p,t,ui,d);

DrawSpline(n,p,P,ui);

[epsilon,e] = shujudianwucha(M,N,P,d);