clear all;
% 实验1
x=-2:0.05:2-0.05;
y=sin(2*x)+2*exp(-30*x.^2)+2;
% % 实验2
% x=0:0.4:16*pi;
% y=sin(x);
% % 实验3
% t=-4*pi:0.2:4*pi;
% x=sin(0.75*t);
% y=sin(t);


d=[x;y]';%给定数据点,di=(xi,yi) i=1~M+1
dy = gradient(y);
dx = gradient(x);
l=[-dy;dx]';%在数据点处的法向量约束条件li i=1~M+1

M=length(x)-1; %论文中是0~M，所以总数length=M+1

t=canshuhua(M,d);
p=3;%B样条次数p=3，控制顶点n个，节点矢量ui i=1~n+p+2==M+1+p-1,
n=M+p;


% ui=jiedianxiangliang(n,p);%方案1 均匀节点向量
% ui=jiedianxiangliang2(n,p,t);%方案2 平均节点向量
ui=jiedianxiangliang3(n,p,t,M);%方案3 皮格尔逼近节点向量


[N,R,P] = kongzhidingdian(M,n,p,t,ui,d);

DrawSpline(n,p,P,ui);
axis([x(1)-1,x(end)+1,min(y)-1,max(y)+1])

[epsilon,e] = shujudianwucha(M,N,P,d);
[sst,es,ev] = faxiangliangwucha(M,n,p,ui,t,l,P);

fprintf('数据点处误差是%d\n法向量误差是%d',e,ev)