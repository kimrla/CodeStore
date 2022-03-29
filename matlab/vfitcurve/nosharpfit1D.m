close all
clear
plan=3;
switch plan
    
    case 1
        
        
        x=0:0.005:1;
        f=1./(0.01+(x-0.3).^2).*(x<0.5)+1./(0.015+(x-0.65).^2).*(x>=0.5);
        x=x';
        f=f';
    case 2
        
        
        x=0:0.02:1;
%         f=100./exp(abs(10*x-5));
%         f=(5*x-2.5).^3+(5*x-2.5).^2+0.5;
        f=((4*x-1).^3/2+(4*x-1).^2/2).*(x<0.5)+((3-4*x).^3/2+(3-4*x).^2/2).*(x>=0.5);
        x=x';
        f=f';
    case 3
        load ECGmoni04401.mat
        x=x'/max(x);
        f=ecg';
        if mod(length(x),2)==1
        x(end)=[];
        f(end)=[];
        end
    case 4
        load ECGmoni02401.mat
        x=x'/max(x);
        f=ecg';
        if mod(length(x),2)==1
        x(end)=[];
        f(end)=[];
        end
    case 5
        load ECG105.mat
        x=linspace(0,1,length(nOP))';
        f=nOP;
    case 6
        load ECG109.mat
        x=linspace(0,1,length(nOP))';
        
        

        f=nOP;
    case 7
        load ECG115.mat
        x=linspace(0,1,length(nOP))';
        f=nOP;
end

% f=f+0.1*rand(size(f));
% f=awgn(f,30,'measured');
% gpoint=[x,f];
% num=length(x);
% tt=linspace(0,1,10*length(f));
% ff=interp1(x,f,tt,'spline');
% for i=2:num-1    
%     [K(i),~] = PJcurvature(gpoint,i);
% end
% for i=2:length(tt)-1
%     KK(i)=PJcurvature([tt',ff'],i);
% end
% K(1)=0;
% K(num)=0;
% phi(1)=pi;
% phi(num)=pi;
% q=zeros(num,1);
% for i=2:num-1
%     r=i+1;
%     l=i-1;
%     
%     if K(i)>K(r) && K(i)>K(l) && K(i)>10*mean(K)
%         %         && phi(i)>pi/6        &&K(i)>3*mean(K)
%         q(i)=1;
%     end
%     
% end
% q(1)=1;
% q(end)=1;
% tzloc=find(q==1);
t=x;
% 
% tzloc=[1;PEAK_1;length(x)];
% tezhengt=t(tzloc);
k=3;
N=5;

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

newt=x;





CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=-1;
Lambda = LSCurFit_V(f,k,N,newt,CList);

curve=LSMatrix_V(k,N,newt)*Lambda;
err=abs(f-curve);
meanerr=mean(err);
maxerr=max(err);

for i=1:size(CList,1)
    for j=0:1
        tempCList=CList;
        tempCList(i,3)=j;
        tempLambda=LSCurFit_V(f,k,N,newt,tempCList);
        tempcurve=LSMatrix_V(k,N,newt)*tempLambda;
        temperr=abs(f-tempcurve);
        tempmeanerr=mean(temperr);
        
        if  tempmeanerr<meanerr
            CList(i,3)=j;
            Lambda=tempLambda;            
            meanerr=tempmeanerr;            
        end
    end
end




%% 本方法最终重构结果
figure
curve=LSMatrix_V(k,N,newt)*Lambda;
err=abs(f-curve);
maxerr=max(err);%最大误差
meanerr=mean(err);%平均误差
stderr=std(err);%标准差
MSE=immse(f,curve);%均方误差
CC=XGD(f,curve);
plot(newt,f,'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
plot(newt,curve,'Color',[0 102 153]/255,'LineWidth',3)
plot(newt,err,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('最终重构结果','fontsize', 15, 'fontname', '微软雅黑')
legend({'原始信号','重构信号','绝对误差'},'fontsize', 15, 'fontname', '微软雅黑')
snr1=SNR(f,curve);

% plot(gpoint(locs,1),gpoint(locs,2),'s','MarkerSize',10)
% plot(gpoint(tlist,1),gpoint(tlist,2),'s','MarkerSize',10)
duandian=LSMatrix_V(k,N,CList(:,1))*Lambda;
plot(CList(:,1),duandian,'s','MarkerSize',10)



%% db4小波变换对比实验
figure
[cdb4,ldb4] = wavedec(f,3,'db4');
rdb4=cdb4;
rdb4(ldb4(1)+1:end) = 0;
iwpdb4 = waverec(rdb4,ldb4,'db4');
errdb4=abs(f-iwpdb4);
maxerrdb4=max(errdb4);%最大误差
meanerrdb4=mean(errdb4);%平均误差
stderrdb4=std(errdb4);%标准差
MSEdb4=immse(f,iwpdb4);%均方误差
CCdb4=XGD(f,iwpdb4);
plot(f,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(iwpdb4,'Color',[0 102 153]/255,'LineWidth',3)
plot(errdb4,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(ldb4(1)),'项db4小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
snrdb4=SNR(f,iwpdb4);
%% db2小波变换对比实验
figure
[cdb2,ldb2] = wavedec(f,3,'db2');
rdb2=cdb2;
rdb2(ldb2(1)+1:end) = 0;
iwpdb2 = waverec(rdb2,ldb2,'db2');
errdb2=abs(f-iwpdb2);
maxerrdb2=max(errdb2);%最大误差
meanerrdb2=mean(errdb2);%平均误差
stderrdb2=std(errdb2);%标准差
MSEdb2=immse(f,iwpdb2);%均方误差
CCdb2=XGD(f,iwpdb2);
plot(f,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(iwpdb2,'Color',[0 102 153]/255,'LineWidth',3)
plot(errdb2,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(ldb2(1)),'项db2小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
snrdb2=SNR(f,iwpdb2);
%% Fourier变换对比试验
figure
nd=64;
z = frdescp1d(f);
P_FFT = ifrdescp1d(z,nd);
errf=abs(f-P_FFT);
maxerrf=max(errf);%最大误差
meanerrf=mean(errf);%平均误差
stderrf=std(errf);%标准差
MSEf=immse(f,P_FFT);%均方误差
CCf=XGD(f,P_FFT);
plot(f,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(P_FFT,'Color',[0 102 153]/255,'LineWidth',3)
plot(errf,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(nd),'项Fourier变换'],'fontsize', 15, 'fontname', '微软雅黑')
snrf=SNR(f,P_FFT);
%% 累积归一化能量图
figure
energyv=cumsum(sort(abs(Lambda(5:64)).^2,'descend'));
normenergyv=energyv/energyv(end);
% energyf=cumsum(sort(abs(z(4:64)).^2,'descend'));
energyf=cumsum(sort(abs(z([length(z)/2-32:length(z)/2-3,length(z)/2+3:length(z)/2+32])).^2,'descend'));
normenergyf=energyf/energyf(end);
energydb4=cumsum(sort(abs(cdb4(5:64)).^2,'descend'));
nomenergydb4=energydb4/energydb4(end);
energydb2=cumsum(sort(abs(cdb2(5:64)).^2,'descend'));
nomenergydb2=energydb2/energydb2(end);
% energystl=[0;cumsum(sort(sortedlam.^2,'descend'))];
% normenergystl=energystl/energystl(end);
hold on
plot(normenergyv,'LineWidth',3)
plot(normenergyf,'LineWidth',3)
plot(nomenergydb4,'LineWidth',3)
plot(nomenergydb2,'LineWidth',3)
% plot(normenergystl(1:256),'LineWidth',1.5)
legend({'V系统正交变换','Fourier变换','db4小波变换','db2小波变换'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
xlim([0 64])
ylim([0 1.2])
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('累积归一化频谱能量图','fontsize', 15, 'fontname', '微软雅黑')
%%
function [phi,K]=qulv(gpoint,i,num)
r=i+1;
l=i-1;

v1=gpoint(i,:)-gpoint(l,:);
v2=gpoint(r,:)-gpoint(i,:);
v3=gpoint(r,:)-gpoint(l,:);
phi=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
K=2*sin(phi)/norm(v3);
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