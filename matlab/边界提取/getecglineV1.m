clear
clc
close all;
%------ 具体数据------------------------------------------------------
PATH= 'C:\CodeStore\matlab\mit-bih-arrhythmia-database-1.0.0'; % path, 这里就是写刚才你保存的数据地址
ecgname='100';
HEADERFILE=[ecgname,'.hea'];      % 文件格式为文本格式
ATRFILE= [ecgname,'.atr'];         % attributes-file 文件以二进制格式
DATAFILE=[ecgname,'.dat'];         % data-file
SAMPLES2READ=3600;         %  读取的数据样本点数为660000
% in case of more than one signal:
% 2*SAMPLES2READ samples are read
%------ LOAD HEADER DATA --------------------------------------------------
signald= fullfile(PATH, DATAFILE);            % data in format 212
fid2=fopen(signald,'r');
A= fread(fid2, [3, SAMPLES2READ], 'uint8')';  % matrix with 3 rows, each 8 bits long, = 2*12bit
fclose(fid2);
%=----------------------载入二进制数据-----------------------------------------
M2H= bitshift(A(:,2), -4);          %字节向右移四位，即取字节的高四位
M1H= bitand(A(:,2), 15);            %取字节的低四位
M( : , 1)= bitshift(M1H,8)+ A(:,1); %低四位向左移八位
M( : , 2)= bitshift(M2H,8)+ A(:,3); %高四位向左移八位
M = M-1024;                               %这个M就是咱们解码出来的数据
nOPnum=1024;%截取数据个数
stp=120;%截取数据点起点
ORIGINP=M(stp:stp+nOPnum-1,1);%取需要预处理的原始数据
plot (ORIGINP)         %绘制一段心电图形
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('原始信号','fontsize', 15, 'fontname', '微软雅黑')
R=100;%中值滤波的窗口W=2*R+1
nOP=myMedfilt(ORIGINP, R);%中值滤波后的结果
% wdenoise%小波去噪
% figure
% plot (XC1) %XC1为去噪后的数据

figure
plot(nOP);
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('中值滤波处理后的信号','fontsize', 15, 'fontname', '微软雅黑')
% [c,l] = wavedec(nOP,3,'coif5');%小波分析滤波 3层分解
% a3 = wrcoef('a',c,l,'coif5');
% figure
% plot(a3);

k=3;
Nmax=9;
num=(k+1)*2^(Nmax-1);%截取数据点位置
tt=linspace(1,2*num-1,num)'/(2*num);
% stp=1109;%截取数据点起点
% P=a3(stp:512+stp-1);
% P=nOP(stp:512+stp-1);
% P=[linspace(1,length(P),length(P))',P];
% t=canshuhua(P);
% t=linspace(0,1,length(nOP))';
% PP=interp1(t,nOP,tt,'linear');%对原始数据插值取点
PP=nOP;
matrixname=['V',num2str(k),'_',num2str(num),'.mat'];
load(matrixname)
N=1;%初始V组数

tol1 = 1e-1; % 允许的最大误差
tol2 = 1e-6; % 要求的最小误差变化
tol3=1e-1;%允许的最大平均误差

err1 = tol1 + 1;%误差初始值
err2 = tol2 + 1;%误差变化初始值
meanerr=tol3+1;%平均误差初始值
meanout=meanerr;
gen=0;%当前迭代数
Gmax=200;%最大迭代数
maxout=tol1+1;%局部片段误差初始值
tpiece=[0,1];%初始参数区间划分
tidx=[1,num];%初始区间划分索引
Lambda1=[];
% figure
% plot(nOP,'.'),hold on;
% plot(PP,'.')
P0=PP;
while(err1 > tol1 && gen<Gmax&&N<Nmax)
    clear fidx
    fidx=find(maxout>tol1);%找出需要加层的区间索引
    % while(meanerr > tol3 && err2 > tol2 && gen<Gmax&&N<Nmax)
    %     clear fidx
    %     fidx=find(meanout>tol3);%找出需要加层的区间索引
    N=N+1;
    gen=gen+1;
    vmat=LSMatrix_V(k,N,tt);
    %     clear CList
    %     CList=zeros(2^(N-1)-1,3);
    %     CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
    %     CList(:,2)=CList(:,1);
    %     CList(:,3)=2;
    %     A = LSMatrix_V(k,N,tt);
    if N==2
        
        temA=O(:,1:(k+1)*2);
        Lambda0=temA'*PP;
    else
        Lambda0=zeros((k+1)*2^(N-2),1);
        temA=O(:,(k+1)*2^(N-2)+1:(k+1)*2^(N-1));
        
        for i=1:length(fidx)
            temP=PP(tt>tpiece(fidx(i))&tt<tpiece(fidx(i)+1),:);
            Lambda0(i:2^(N-2):i+k*2^(N-2),:)=temA(tt>tpiece(fidx(i))&tt<tpiece(fidx(i)+1),i:2^(N-2):i+k*2^(N-2))'*temP;
        end
    end
    Lambda1=[Lambda1;Lambda0];
    %         figure,
    %         plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
    %         VCompose(Lambda1,k,N)
    %         Lambda1 = LSCurFit_L2LC(Lambda1,k,N,CList);
    
    
    
    %         Lambda3=O(:,1:(k+1)*2^(N-1))'*P0;
%             Lambda2=vmat\P;
    PC=O(:,1:(k+1)*2^(N-1))*Lambda1;
    deltaP=P0- PC;%拟合曲线和原始数据点的差向量
    delta = abs(deltaP);%每个点的误差
    meanerr=mean(delta);
    PP=deltaP;
    tmp=max(delta);
    err2=abs(err1-tmp);
    err1=tmp;
    tpiece=linspace(0,1,2^(N-1)+1);%参数分割区间
    tidx=[];
    for i=1:2^(N-1)
        tidx(end+1)=find(tt>=tpiece(i),1);
    end
    i=i+1;
    tidx(end+1)=num;
    
    m=max(diff(tidx));
    clear out
    out(1:m,1:length(tidx)-1)=nan;
    
    for i=1:length(tidx)-2
        out(1:tidx(i+1)-tidx(i),i)=delta(tidx(i):tidx(i+1)-1);
    end
    i=i+1;
    out(1:tidx(i+1)-tidx(i)+1,i)=delta(tidx(i):tidx(i+1));
    maxout=max(out,[],'omitnan');
    meanout=mean(out,'omitnan');
    % %
    figure,
    plot(nOP,'Color',[255 102 102]/255,'MarkerSize',15);hold on
    %     VCompose(Lambda1,k,N)
    
        plot(PC,'Color',[0 102 153]/255,'LineWidth',1.5)
        set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
        title(['第',num2str(N),'代变换结果'],'fontsize', 15, 'fontname', '微软雅黑')
%     LSVCompose(Lambda1,k,Nmax,O)
    
%         figure,
%         plot(P0,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
%         LSVCompose(Lambda2,k,Nmax,O)
    %     VCompose(Lambdatrapz,k,N)
end

% CList=zeros(2^(N-1)-1,3);
% CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
% CList(:,2)=CList(:,1);
% CList(:,3)=2;
% Lambda1 = LSCurFit_L2LC(Lambda1,k,N,CList);
% figure,
% plot(P0,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
% LSVCompose(Lambda1,k,Nmax,O)
%% 本方法最终重构结果
figure
curve=O(:,1:size(Lambda1,1))*Lambda1;
err=abs(nOP-curve);
maxerr=max(err);%平均误差
meanerr=mean(err);%最大误差
stderr=std(err);%标准差
MSE=immse(nOP,curve);%均方误差
plot(nOP,'Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(curve,'Color',[0 102 153]/255,'LineWidth',1.5)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('最终重构结果','fontsize', 15, 'fontname', '微软雅黑')
%% 本方法将基函数系数排列取有限个最大值重构结果
sortedlam=sort(Lambda1);






%% 本方法取低频基函数重构结果
figure
cutnum=256;
curve1=O(:,1:cutnum)*Lambda1(1:cutnum,:);
errs=abs(nOP-curve1);
maxerrs=max(errs);%平均误差
meanerrs=mean(errs);%最大误差
stderrs=std(errs);%标准差
MSEs=immse(nOP,curve1);%均方误差
plot(nOP,'Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(curve1,'Color',[0 102 153]/255,'LineWidth',1.5)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title(['前',num2str(cutnum),'项V基函数重构结果'],'fontsize', 15, 'fontname', '微软雅黑')

%% db4小波变换对比实验
figure
[cdb4,ldb4] = wavedec(nOP,4,'db4');
cdb4(ldb4(4)+1:end) = 0;
iwpdb4 = waverec(cdb4,ldb4,'db4');
errdb4=abs(nOP-iwpdb4);
maxerrdb4=max(errdb4);%平均误差
meanerrdb4=mean(errdb4);%最大误差
stderrdb4=std(errdb4);%标准差
MSEdb4=immse(nOP,iwpdb4);%均方误差
plot(nOP,'Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(iwpdb4,'Color',[0 102 153]/255,'LineWidth',1.5)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(ldb4(4)),'项db4小波变换'],'fontsize', 15, 'fontname', '微软雅黑')

%% db2小波变换对比实验
figure
[cdb2,ldb2] = wavedec(nOP,4,'db2');
cdb2(ldb2(4)+1:end) = 0;
iwpdb2 = waverec(cdb2,ldb2,'db2');
errdb2=abs(nOP-iwpdb2);
maxerrdb2=max(errdb2);%平均误差
meanerrdb2=mean(errdb2);%最大误差
stderrdb2=std(errdb2);%标准差
MSEdb2=immse(nOP,iwpdb2);%均方误差
plot(nOP,'Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(iwpdb2,'Color',[0 102 153]/255,'LineWidth',1.5)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(ldb2(4)),'项db2小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
%% Fourier变换对比试验
figure
nd=256;
z = frdescp1d(nOP);
P_FFT = ifrdescp1d(z,nd);
errf=abs(nOP-P_FFT);
maxerrf=max(errf);%平均误差
meanerrf=mean(errf);%最大误差
stderrf=std(errf);%标准差
MSEf=immse(nOP,P_FFT);%均方误差
plot(nOP,'Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(P_FFT,'Color',[0 102 153]/255,'LineWidth',1.5)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(nd),'项傅里叶变换'],'fontsize', 15, 'fontname', '微软雅黑')