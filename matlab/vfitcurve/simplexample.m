%%4.5论文2的简单例子
close all
clear
num=201;
x=linspace(0,1,num)';
f=(0.005./(0.01+(x-0.25).^2)-0.005/(0.01+0.25^2)).*(x<=0.5)+(1./exp(abs(20*x-15))-1/exp(5)).*(x>0.5);
for i=2:num-1
    k(i)=PJcurvature([x,f],i);
end
[~,tzloc]=findpeaks(k,'MinPeakHeight',mean(k));
tzloc=[1 tzloc num];
tzx=x(tzloc);

k=3;
N=5;

CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=2;
Lambda = LSCurFit_V(f,k,N,x,CList);

curve=LSMatrix_V(k,N,x)*Lambda;
err=vecnorm(f-curve,2,2);
meanerr=mean(err);
maxerr=max(err);
plot(f,'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
plot(curve,'Color',[0 102 153]/255,'LineWidth',3)

for i=2:length(tzloc)-1
    for j=0:1
        tempCList=CList;
        tempCList(ismember(CList(:,1),tzx(i)),3)=j;
        tempLambda=LSCurFit_V(f,k,N,x,tempCList);
        tempcurve=LSMatrix_V(k,N,x)*tempLambda;
        temperr=vecnorm(f-tempcurve,2,2);
        tempmeanerr=mean(temperr);
        
        if  tempmeanerr<meanerr
            CList(ismember(CList(:,1),tzx(i)),3)=j;
            Lambda=tempLambda;            
            meanerr=tempmeanerr;            
        end
    end
end
figure
curve=LSMatrix_V(k,N,x)*Lambda;
err=vecnorm(f-curve,2,2);
maxerr=max(err);%最大误差
meanerr=mean(err);%平均误差
stderr=std(err);%标准差
MSE=immse(f,curve);%均方误差
plot(f,'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
plot(curve,'Color',[0 102 153]/255,'LineWidth',3)






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