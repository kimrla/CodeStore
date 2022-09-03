%%4.5论文2的简单例子
close all
clear
example=1;
switch example
    case 1
num=101;
x=linspace(0,1,num)';
% f=(0.004./(0.01+(x-0.25).^2)-0.004/(0.01+0.25^2)).*(x<=0.5)+(1./exp(abs(20*x-15))-1/exp(5)).*(x>0.5);
% f=(-8*x.^2+4*x).*(x<0.25)+(64*x.^3-72*x.^2+24*x-2).*(x>=0.25 & x<0.5) ...
%     +((4*x-2).^3/2+(4*x-2).^2/2).*(x>=0.5 & x<0.75)+((4-4*x).^3/2+(4-4*x).^2/2).*(x>=0.75);
f=(6*x.^3-9*x.^2+3*x).*(x<0.5) ...
    +((4*x-2).^3/2+(4*x-2).^2/2).*(x>=0.5 & x<0.75)+((4-4*x).^3/2+(4-4*x).^2/2).*(x>=0.75);
case 2
        % 实验2 尖锐特征
        x=0:0.002:1-0.002;
        x=x';
        num=length(x);
        f=100./exp(abs(10*x-6))+(10*x-6).^5/500;
end
% load simplexample101.mat
% f=curve;
for i=2:num-1
    kappa(i)=PJcurvature([x,f],i);
end
P=[x,f];
[~,tzloc]=findpeaks(kappa,'MinPeakHeight',mean(kappa));
tzloc=[1 tzloc num];
t=x;
tezhengt=t(tzloc);

k=3;
N=4;

for i=1:length(tezhengt)
    newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),N-1));
end
while length(newtezhengt)-length(unique(newtezhengt))
    N=N+1;
    for i=1:length(tezhengt)
        newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),N-1));
    end
end

newt(tzloc)=newtezhengt;
for i=2:length(tzloc)
    j=tzloc(i-1)+1:tzloc(i)-1;
    newt(j)=(t(j)-t(tzloc(i-1)))*(newtezhengt(i)-newtezhengt(i-1))/(tezhengt(i)-tezhengt(i-1))+newt(tzloc(i-1));
end

newt=newt';
% newt=t;
CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=2;

CList(ismember(CList(:,1),newtezhengt(2)),3)=2;
CList(ismember(CList(:,1),newtezhengt(3)),3)=0;
CList(ismember(CList(:,1),newtezhengt(4)),3)=0;
Lambda = LSCurFit_V(P,k,N,newt,CList);
% Lambda=LSMatrix_V(k,N,newt)\f;
curve=LSMatrix_V(k,N,newt)*Lambda;
err=abs(f-curve);
meanerr=mean(err);
maxerr=max(err);
plot(x,f,'.','Color',[255 102 102]/255,'MarkerSize',10,'linewidth',1.5);hold on
plot(x(tzloc),f(tzloc),'s','MarkerSize',10,'MarkerFaceColor',[0 102 153]/255)
% VCompose(Lambda,k,N);
plot(x,curve(:,2),'Color',[0 102 153]/255,'LineWidth',3)
legend({'原始数据','特征点','拟合曲线'},'fontsize', 15, 'fontname', '微软雅黑')

set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0])
for i=2:length(tzloc)-1
    for j=0
        tempCList=CList;
        tempCList(ismember(CList(:,1),newtezhengt(i)),3)=j;
        tempLambda=LSCurFit_V(f,k,N,newt,tempCList);
        tempcurve=LSMatrix_V(k,N,newt)*tempLambda;
        temperr=abs(f-tempcurve);
        tempmeanerr=mean(temperr);
        
        if  tempmeanerr<meanerr
            CList(ismember(CList(:,1),newtezhengt(i)),3)=j;
            Lambda=tempLambda;            
            meanerr=tempmeanerr;            
        end 
    end    
end
figure
curve=LSMatrix_V(k,N,newt)*Lambda;
err=vecnorm(f-curve,2,2);
maxerr=max(err);%最大误差
meanerr=mean(err);%平均误差
stderr=std(err);%标准差
MSE=immse(f,curve);%均方误差
plot(f,'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
plot(curve,'Color',[0 102 153]/255,'LineWidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0],'Ylim',[-0.2 1.2])





function [kappa,norm_k] = PJcurvature(gpoint,i)
    x = gpoint(i-1:i+1,1);
    y = gpoint(i-1:i+1,2);
    t_a = norm([x(2)-x(1),y(2)-y(1)]);
    t_b = norm([x(3)-x(2),y(3)-y(2)]);
    
    M =[[1, -t_a, t_a^2];
        [1, 0,    0    ];
        [1,  t_b, t_b^2]];

    a = M\x;
    b = M\y;

    kappa  = 2.*abs(a(3)*b(2)-b(3)*a(2)) / (a(2)^2.+b(2)^2.)^(1.5);
    norm_k =  [b(2),-a(2)]/sqrt(a(2)^2.+b(2)^2.);
end
function y=xiaoshu2dec(x,N)
y=zeros(1,N+1);
tempx=x;
for i=1:N+1
    
    y(i)=floor(tempx*2);
    tempx=tempx*2-y(i);
end
if y(end)==1
    y(N)=y(N)+1;
    y(end)=[];
end
for i=fliplr(2:N)
    if y(i)==2
        y(i)=0;
        y(i-1)=y(i-1)+1;
    end
end
% if y(1)==2
%     for i=1:N
%         y(i)=floor(x*2);
%         x=x*2-y(i);
%     end
% end
end

function y=dec2xiaoshu(x)
for i=1:length(x)
    temp(i)=x(i)*2^(-i);
end
y=sum(temp);
end