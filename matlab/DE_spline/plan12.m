clear all;
% 实验1
x=0:0.01:1-0.01;
f=90./(1+exp(-100*(x-0.4)));
n=8;%控制顶点个数-1
N=length(x);%采样点个数

% % 实验2
% x=0:0.05:10;
% f=100./exp(abs(x-5))+(x-5).^2/500;
% n=13;
% N=length(x);%采样点个数

f_=f+normrnd(0,0.01,1,N);



plot(x,f_,"*")
hold on
% plot(x,f)
d=[x;f_]';%给定数据点,di=(xi,yi) i=1~M+1




M=length(x)-1; %论文中是0~M，所以总数length=M+1

t=canshuhua(M,d);
p=3;%B样条次数p=3，控制顶点n+1个，节点矢量ui i=1~n+p+2,




ui=jiedianxiangliang(n,p);%方案1 均匀节点向量


% u=0:0.01:1;
% for j=1:n+1
%     for i=1:length(u)
%         Njp_u(i,j) = Njp(j, p , u(i), ui);  
%     end
%     subplot(5,5,j),plot(u,Njp_u(:,j));
% end


[N,R,P] = kongzhidingdian(M,n,p,x,ui,d);

DrawSpline(n,p,P,ui);

[epsilon,e] = shujudianwucha(M,N,P,d);