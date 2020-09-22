clear all;
% % 实验3
% [x,y]=meshgrid(0:0.01:1);
% f=74*(x>2/3&x<=1&y>=0&y<=1/3)+76*((x>=1/3&x<=2/3&y>=0&y<=1/3)|(x>2/3&x<=1&y>1/3&y<=2/3))+...
%     78*((x>=0&x<=1/3&y>=0&y<=1/3)|(x>1/3&x<=2/3&y>1/3&y<=2/3)|(x>2/3&x<=1&y>2/3&y<=1))+...
%     80*((x>=0&x<=1/3&y>1/3&y<=2/3)|(x>1/3&x<=2/3&y>2/3&y<=1))+...
%     82*(x>=0&x<=1/3&y>2/3&y<=1);
% nu=15;
% nv=15;

% 实验4
[x,y]=meshgrid(0:0.01:1);
f=exp(-0.04*((80.*x-40).^2+(90.*y-45).^2).^0.5).*cos(0.15*((80.*x-40).^2+(90.*y-45).^2).^0.5);
nu=9;
nv=9;

plot3(x,y,f,"*")
mesh(x,y,f)
hold on
dx=[x;f]';%给定数据点,di=(xi,yi) i=1~M+1
dy=[y;f]';

Mx=length(x)-1; %论文中是0~M，所以总数length=M+1
My=length(y)-1;
tx=canshuhua(Mx,dx);
ty=canshuhua(My,dy);
p=3;%B样条次数p=3，控制顶点n+1个，节点矢量ui i=1~n+p+2,




ui=jiedianxiangliang(nu,p);% 均匀节点向量
vi=jiedianxiangliang(nv,p);

[Nu,~,Pu] = kongzhidingdian(Mx,nu,p,tx,ui,dx);
[Nv,~,Pv] = kongzhidingdian(My,nv,p,ty,vi,dy);
DrawSpline(nu,p,Pu,ui);
DrawSpline(nv,p,Pv,vi);
% [epsilon,e] = shujudianwucha(M,N,P,d);