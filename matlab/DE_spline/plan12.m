clear all;
% % 实验1
% x=0:0.025:1;
% y=90./(1+exp(-100*(x-0.4)));
% n=8;

% 实验2
x=0:0.05:10;
y=100./exp(abs(x-5))+(x-5).^2/500;
n=13;

plot(x,y,"*")
hold on
d=[x;y]';%给定数据点,di=(xi,yi) i=1~M+1


M=length(x)-1; %论文中是0~M，所以总数length=M+1

t=canshuhua(M,d);
p=3;%B样条次数p=3，控制顶点n+1个，节点矢量ui i=1~n+p+2,




ui=jiedianxiangliang(n,p);%方案1 均匀节点向量


[N,R,P] = kongzhidingdian(M,n,p,t,ui,d);

DrawSpline(n,p,P,ui);

[epsilon,e] = shujudianwucha(M,N,P,d);