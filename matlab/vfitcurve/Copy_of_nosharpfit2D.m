% 不提取特征，直接通过加密采样和堆叠基函数的算法
close
clear
plan=2;
switch plan    
    case 1
        load bird200.mat
        N=4;
        rd=3;
    case 2
        load fire500.mat
        N=10;
        rd=2;
        gpoint(end+1,:)=gpoint(1,:);
    case 3
        load yezi600.mat
        N=4;
        rd=1;
        gpoint(end+1,:)=gpoint(1,:);
    case 4
        load shizi1500.mat
        N=4;
        rd=0.1;
    case 5
        load fenghuang2000.mat
        N=10;
%         N=3;
        rd=0.01;
    case 6
        load hudie3000.mat
%         N=10;
        N=3;
        rd=0.01;
    case 7
        load niao21000.mat
%         N=7;
        N=3;
        rd=0.2;
    case 8
        load huacao4-1500.mat
        
%         N=10;
        N=3;
        rd=0.1;
    case 9
        load G-200.mat
%         N=6;
        N=2;
        rd=1;
        gpoint(end+1,:)=gpoint(1,:);
end
P=gpoint;
num=length(P);
k=3;
t=canshuhua(P);
for i=2:num-1
%     [phi(i),K(i)]=qulv(gpoint,i,num);
    [K(i),~] = PJcurvature(gpoint,i);
end
K(num)=K(1);
q=zeros(num,1);
for i=1:num
    r=i+1;
    l=i-1;
    if i==1
        l=num-1;
    elseif i==num
        r=2;
    end
    if K(i)>K(r) && K(i)>K(l) && K(i)>mean(K)
%         phi(i)>pi/6 && (phi(i)>pi/6)*(phi(r)>pi/6)*(phi(l)>pi/6)==0
        %         && phi(i)>pi/6        &&K(i)>3*mean(K)
        q(i)=1;
    end
    
end
q(1)=1;
q(num)=1;
tzloc=find(q==1);
% tezhengt=t(tzloc);
% 
% for i=1:length(tezhengt)
%     newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),N-1));
% end
% while length(newtezhengt)-length(unique(newtezhengt))
%     N=N+1;
%     for i=1:length(tezhengt)
%         newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),N-1));
%     end
% end
% 
% newt(tzloc)=newtezhengt;
% for i=2:length(tzloc)
%     j=tzloc(i-1)+1:tzloc(i)-1;
%     newt(j)=(t(j)-t(tzloc(i-1)))*(newtezhengt(i)-newtezhengt(i-1))/(tezhengt(i)-tezhengt(i-1))+newt(tzloc(i-1));
% end
% 
% newt=newt';




if (k+1)*2^(N-1)>length(P)
    tt=unique([linspace(0,1,(k+1)*2^(N-1))';t]);
    P=interp1(t,P,tt,'linear');
    t=tt;
end
newt=t;
newtezhengt=t(tzloc);
% f=f+0.1*rand(size(f));
CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=2;
Lambda = LSCurFit_V(P,k,N,newt,CList);

curve=LSMatrix_V(k,N,newt)*Lambda;
err=vecnorm(P-curve,2,2);
meanerr=mean(err);
maxerr=max(err);

for i=2:length(tzloc)-1
    for j=0:1
        tempCList=CList;
        tempCList(ismember(CList(:,1),newtezhengt(i)),3)=j;
        tempLambda=LSCurFit_V(P,k,N,newt,tempCList);
        tempcurve=LSMatrix_V(k,N,newt)*tempLambda;
        temperr=vecnorm(P-tempcurve,2,2);
        tempmeanerr=mean(temperr);
        
        if  tempmeanerr<meanerr
            CList(ismember(CList(:,1),newtezhengt(i)),3)=j;
            Lambda=tempLambda;            
            meanerr=tempmeanerr;            
        end
    end
end
LambdaS=Lambda;
LambdaS(all(abs(Lambda)<=5*10^(-2),2),:)=0;
figure
curveS=LSMatrix_V(k,N,newt)*LambdaS;
errS=vecnorm(P-curveS,2,2);
maxerrS=max(errS);%最大误差
meanerrS=mean(errS);%平均误差
stderrS=std(errS);%标准差
MSES=immse(P,curveS);%均方误差
plot(gpoint(:,1),gpoint(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
plot(curveS(:,1),curveS(:,2),'Color',[0 102 153]/255,'LineWidth',3)
LambdaS(all(LambdaS~=0,2),3)=find(all(abs(LambdaS)~=0,2));
LambdaS(all(LambdaS==0,2),:)=[];


%% 本方法最终重构结果
figure
curve=LSMatrix_V(k,N,newt)*Lambda;
err=vecnorm(P-curve,2,2);
maxerr=max(err);%最大误差
meanerr=mean(err);%平均误差
stderr=std(err);%标准差
MSE=immse(P,curve);%均方误差

plot(gpoint(:,1),gpoint(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
plot(curve(:,1),curve(:,2),'Color',[0 102 153]/255,'LineWidth',3)
% VCompose(Lambda,k,N)
% plot(t,err,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('最终重构结果','fontsize', 15, 'fontname', '微软雅黑')
legend({'原始信号','重构信号'},'fontsize', 15, 'fontname', '微软雅黑')


% plot(P(tzloc,1),P(tzloc,2),'s','MarkerSize',10)
% plot(P(tzloc,1),P(tzloc,2),'s','MarkerSize',10)
duandian=LSMatrix_V(k,N,CList(:,1))*Lambda;
% plot(duandian(:,1),duandian(:,2),'s','MarkerSize',10)

function [phi,K]=qulv(gpoint,i,num)
r=i+1;
l=i-1;
if i==1
    l=num-1;
elseif i==num
    r=2;
end
v1=gpoint(i,:)-gpoint(l,:);
v2=gpoint(r,:)-gpoint(i,:);
v3=gpoint(r,:)-gpoint(l,:);
phi=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
K=2*sin(phi)/norm(v3);
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