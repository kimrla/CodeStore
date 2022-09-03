%% 对ECG局部增量式信号变换
clear
clc
close all;
%------ 具体数据------------------------------------------------------
PATH= 'C:\CodeStore\matlab\mit-bih-arrhythmia-database-1.0.0'; % path, 这里就是写刚才你保存的数据地址
ecgname='105'; %105 198+1024  %109 775+1024 %115 958+1024 %201 1049+1024不好 %100 1049+1024 不好 %104 1049+1024不好
nOPnum=1024;%截取数据个数
stp=198;%截取数据点起点
HEADERFILE=[ecgname,'.hea'];      % 文件格式为文本格式
ATRFILE= [ecgname,'.atr'];         % attributes-file 文件以二进制格式
DATAFILE=[ecgname,'.dat'];         % data-file
%% 读取头文件
% 通过函数 fullfile 获得头文件的完整路径
signalh= fullfile(PATH, HEADERFILE);
% 打开头文件，其标识符为 fid1 ，属性为'r'--“只读”
fid1=fopen(signalh,'r');
% 读取头文件的第一行数据，字符串格式
z= fgetl(fid1);
% 按照格式 '%*s %d %d %d' 转换数据并存入矩阵 A 中
A= sscanf(z, '%*s %d %d %d',[1,3]);

nosig= A(1);    % 信号通道数目
sfreq=A(2);     % 数据采样频率
SAMPLES2READ = 10*sfreq;    %取十秒数据

for k=1:nosig           % 读取每个通道信号的数据信息
    z= fgetl(fid1);
    A= sscanf(z, '%*s %d %d %d %d %d',[1,5]);
    dformat(k)= A(1);       % 信号格式; 这里只允许为 212 格式
    gain(k)= A(2);          % 每 mV 包含的整数个数
    bitres(k)= A(3);        % 采样精度（位分辨率）
    zerovalue(k)= A(4);     % ECG 信号零点相应的整数值
    firstvalue(k)= A(5);    % 信号的第一个整数值 (用于偏差测试)
end
fclose(fid1);
clear A;
%% 读取信号数据
%.dat文件的数据格式读取为每行三个字节,即三个八位的二进制数字,其内容含义为
%      0000 0000  ||             0000 0000              ||  0000 0000
%sign1(L)低八位信息||左四位sign2(R)信息,右四位sign1(L)信息||sign2(R)低八位信息
%将第二字节的信息处理后,后四位移至第一字节最左位即得到完整的sign1
%将第二字节的信息处理后,前四位移至第一字节最左位即得到完整的sign2.
signald = fullfile(PATH , DATAFILE);
fid2 = fopen(signald,'r');
A= fread(fid2, [3, SAMPLES2READ], 'uint8')';
fclose(fid2);

%对第二字节做左位移运算，位移距离-4
%得到第二字节左四位，即sign2的高四位，包括符号位,右高信息
M_R_H = bitshift(A(:,2), -4);
%对第二字节和1111做与运算，
%保留第二字节右四位，即sign1的低四位，包括符号位,左高信息
M_L_H = bitand(A(:,2), 15);
%对第二字节和1000做与运算，
%保留第二字节右边第四位，获取sign2符号位,并向左位移九位，与整体sign1进行运算
PRL=bitshift(bitand(A(:,2),8),9);
%对第二字节和10000000做与运算，
%保留第二字节右边第四位，获取sign1符号位，并向左位移5位，与整体sign2进行运算
PRR=bitshift(bitand(A(:,2),128),5);

%M矩阵为sign1，2的存储矩阵，存储100.dat处理后数据
%将sign1(L)高位移至sign1低位前(A(:,1))
%将sign2(R)高位移至sign2低位前(A(:,3))
%最后将信号符号位信息去掉
M( : , 1)= bitshift(M_L_H,8)+ A(:,1)-PRL;
M( : , 2)= bitshift(M_R_H,8)+ A(:,3)-PRR;

%将sign的数值与零点做减法得到正负值
%再将得到的具有正负性的值与每mV的整数值相除，即得到电压多少mV
M( : , 1)= (M( : , 1)- zerovalue(1))/gain(1);
M( : , 2)= (M( : , 2)- zerovalue(2))/gain(2);
%将我们设定的采样个数除以频率即得到这段样品的测定时间。
TIME =(0:(SAMPLES2READ-1))/sfreq;
%释放变量
clear A M_R_H M_L_H PRR PRL;
%% 读取注释文件
atrd= fullfile(PATH, ATRFILE);      % attribute file with annotation data
fid3=fopen(atrd,'r');
A= fread(fid3, [2, inf], 'uint8')';
fclose(fid3);
ATRTIME=[];
ANNOT=[];
sa=size(A);
saa=sa(1);
i=1;
while i<=saa
    annoth=bitshift(A(i,2),-2);
    if annoth==59
        ANNOT=[ANNOT;bitshift(A(i+3,2),-2)];
        ATRTIME=[ATRTIME;A(i+2,1)+bitshift(A(i+2,2),8)+...
                bitshift(A(i+1,1),16)+bitshift(A(i+1,2),24)];
        i=i+3;
    elseif annoth==60   % nothing to do!
    elseif annoth==61   % nothing to do!
    elseif annoth==62   % nothing to do!
    elseif annoth==63
        hilfe=bitshift(bitand(A(i,2),3),8)+A(i,1);
        hilfe=hilfe+mod(hilfe,2);
        i=i+hilfe/2;
    else
        ATRTIME=[ATRTIME;bitshift(bitand(A(i,2),3),8)+A(i,1)];
        ANNOT=[ANNOT;bitshift(A(i,2),-2)];
    end
   i=i+1;
end
%ANNOT存储了该记录中各心拍的人工标注类型
ANNOT(length(ANNOT))=[];       % last line = EOF (=0)
%ATRTIME文件存储了各心拍的人工标注位置
ATRTIME(length(ATRTIME))=[];   % last line = EOF
clear A;
ATRPOSITION= (cumsum(ATRTIME));
ATRTIME= (cumsum(ATRTIME))/sfreq;

ind= find(ATRTIME <= TIME(end));
ATRTIMED= ATRTIME(ind);
ANNOT=round(ANNOT);
ANNOTD= ANNOT(ind);
%% 处理信号数据
% nOPnum=1024;%截取数据个数
% stp=198;%截取数据点起点
ORIGINP=M(stp:stp+nOPnum-1,1);%取需要预处理的原始数据
plot (ORIGINP)         %绘制一段心电图形
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('原始信号','fontsize', 15, 'fontname', '微软雅黑')

nOP=ORIGINP;

% figure %中值滤波
% % R=100;%中值滤波的窗口W=2*R+1
% % nOP=myMedfilt(ORIGINP, R);%中值滤波后的结果
% plot(nOP);
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% title('中值滤波处理后的信号','fontsize', 15, 'fontname', '微软雅黑')


figure
% [c,l] = wavedec(nOP,2,'coif5');%小波分析滤波 3层分解
% a3 = wrcoef('a',c,l,'coif5');
% plot(a3);
% nOP=a3;

% [XC1,CXC1,LXC1,PERF0_1,PERF2_1]=wdenoise(nOP);%小波去噪
% nOP=XC1';
% plot (XC1) %XC1为去噪后的数据
PEAK_1=ATRPOSITION(ATRPOSITION>=stp & ATRPOSITION<stp+nOPnum)-stp+1;
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('小波滤波处理后的信号','fontsize', 15, 'fontname', '微软雅黑')

%% 简单例子

close all
clear all
% x=linspace(1,2*128-1,128)'/(2*128)';
x=linspace(0,1,256)';
f=3.*(x<0.25)+2*sin(30*pi*x).*(0.25<=x & x<0.5)+1./exp(abs(20*x-15)).*(x>=0.5);
nOP=f;
plot(x,nOP,'-.','Color',[255 102 102]/255,'MarkerSize',5,'LineWidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0])
%% 局部增量式正交V变换
tic
k=3;
Nmax=7;
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
tol3=1e-3;%允许的最大平均误差

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
temA=O(:,1:k+1);
Lambda0=temA'*PP;
PC=O(:,1:(k+1)*2^(N-1))*Lambda0;
deltaP=P0- PC;%拟合曲线和原始数据点的差向量
delta = abs(deltaP);%每个点的误差
figure,
plot(x,nOP,'-.','Color',[255 102 102]/255,'MarkerSize',15,'LineWidth',3);hold on
plot(x,PC,'Color',[0 102 153]/255,'LineWidth',3)
plot(x,delta,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0])
%     legend({'原始信号','重构信号','绝对误差'},'fontsize', 15, 'fontname', '微软雅黑')
legend({'原始数据','拟合曲线','绝对误差'},'location','best','fontsize', 15, 'fontname', '微软雅黑')


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
            Lambda0(fidx(i):2^(N-2):fidx(i)+k*2^(N-2),:)=temA(tt>tpiece(fidx(i))&tt<tpiece(fidx(i)+1),fidx(i):2^(N-2):fidx(i)+k*2^(N-2))'*temP;
        end
    end
    Lambda1=[Lambda1;Lambda0];
    
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
    
    figure,
    plot(x,nOP,'-.','Color',[255 102 102]/255,'MarkerSize',15,'LineWidth',3);hold on
    plot(x,PC,'Color',[0 102 153]/255,'LineWidth',3)
    plot(x,delta,'linewidth',3)
    set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0])
%     legend({'原始信号','重构信号','绝对误差'},'fontsize', 15, 'fontname', '微软雅黑')
    legend({'原始数据','拟合曲线','绝对误差'},'fontsize', 15, 'fontname', '微软雅黑')
%     title(['第',num2str(N),'代变换结果'],'fontsize', 15, 'fontname', '微软雅黑')
    
end
not0row=sum(any(Lambda1,2));
elapsedTimeV=toc;
drawheat(k,N,Lambda1)
draw2tree(k,N,Lambda1)
drawsumheat(k,N,Lambda1)
%% 本方法最终重构结果
figure
curve=O(:,1:size(Lambda1,1))*Lambda1;
err=abs(nOP-curve);
maxerr=max(err);%最大误差
meanerr=mean(err);%平均误差
stderr=std(err);%标准差
MSE=immse(nOP,curve);%均方误差
CC=XGD(nOP,curve);
plot(nOP,'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
plot(curve,'Color',[0 102 153]/255,'LineWidth',3)
plot(err,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('最终重构结果','fontsize', 15, 'fontname', '微软雅黑')
legend({'原始信号','重构信号','绝对误差'},'fontsize', 15, 'fontname', '微软雅黑')
snr1=SNR(nOP,curve);
%% 本方法将基函数系数排列取有限个最大值重构结果
figure
[~,slidx]=sort(abs(Lambda1),'descend');
firstnum=256;
sortedlam=zeros(length(Lambda1),1);
sortedlam(slidx(1:firstnum))=Lambda1(slidx(1:firstnum));
curvestl=O(:,1:size(Lambda1,1))*sortedlam;
errstl=abs(nOP-curvestl);
maxerrstl=max(errstl);
meanerrstl=mean(errstl);
stderrstl=std(errstl);
MSEstl=immse(nOP,curvestl);
CCstl=XGD(nOP,curvestl);
plot(nOP,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(curvestl,'Color',[0 102 153]/255,'LineWidth',3)
plot(errstl,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title(['系数最大',num2str(firstnum),'项V基函数重构结果'],'fontsize', 15, 'fontname', '微软雅黑')
snrstl=SNR(nOP,curvestl);

%% 本方法取低频基函数重构结果
figure
cutnum=256;
curves=O(:,1:cutnum)*Lambda1(1:cutnum,:);
errs=abs(nOP-curves);
maxerrs=max(errs);%最大误差
meanerrs=mean(errs);%平均误差
stderrs=std(errs);%标准差
MSEs=immse(nOP,curves);%均方误差
CCs=XGD(nOP,curves);
plot(nOP,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(curves,'Color',[0 102 153]/255,'LineWidth',3)
plot(errs,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title(['前',num2str(cutnum),'项V基函数重构结果'],'fontsize', 15, 'fontname', '微软雅黑')
snrs=SNR(nOP,curves);
% %% 连续V系统拟合效果
% figure
% lambdaV=LSMatrix_V(k,7,tt)\nOP;
% curveV=LSMatrix_V(k,7,tt)*lambdaV;
% plot(nOP,'Color',[255 102 102]/255,'MarkerSize',15);hold on
% plot(curveV,'Color',[0 102 153]/255,'LineWidth',1.5)
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% title(['前',num2str(7),'组连续V基函数重构结果'],'fontsize', 15, 'fontname', '微软雅黑')
% 

%% db4小波变换对比实验
tic
figure
[cdb4,ldb4] = wavedec(nOP,2,'db4');
rdb4=cdb4;
rdb4(ldb4(1)+1:end) = 0;
iwpdb4 = waverec(rdb4,ldb4,'db4');
errdb4=abs(nOP-iwpdb4);
maxerrdb4=max(errdb4);%最大误差
meanerrdb4=mean(errdb4);%平均误差
stderrdb4=std(errdb4);%标准差
MSEdb4=immse(nOP,iwpdb4);%均方误差
CCdb4=XGD(nOP,iwpdb4);
plot(nOP,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(iwpdb4,'Color',[0 102 153]/255,'LineWidth',3)
plot(errdb4,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(ldb4(1)),'项db4小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
snrdb4=SNR(nOP,iwpdb4);
elapsedTimedb4=toc;
%% db2小波变换对比实验
figure
[cdb2,ldb2] = wavedec(nOP,2,'db2');
rdb2=cdb2;
rdb2(ldb2(1)+1:end) = 0;
iwpdb2 = waverec(rdb2,ldb2,'db2');
errdb2=abs(nOP-iwpdb2);
maxerrdb2=max(errdb2);%最大误差
meanerrdb2=mean(errdb2);%平均误差
stderrdb2=std(errdb2);%标准差
MSEdb2=immse(nOP,iwpdb2);%均方误差
CCdb2=XGD(nOP,iwpdb2);
plot(nOP,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(iwpdb2,'Color',[0 102 153]/255,'LineWidth',3)
plot(errdb2,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(ldb2(1)),'项db2小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
snrdb2=SNR(nOP,iwpdb2);
%% coif

figure
tic
% nexttile
levelcoif3=3;
cutcoif3level=1;
[ccoif3,lcoif3] = wavedec(P(:,1),levelcoif3,'coif3');
rcoif3=ccoif3;
rcoif3(sum(lcoif3(1:1+cutcoif3level))+1:end) = 0;
iwcoif3 = waverec(rcoif3,lcoif3,'coif3');

elapsedTimecoif3=toc;
energycoif3=iwcoif3.^2;
errcoif3=vecnorm(P-iwcoif3);
maxerrcoif3=max(errcoif3);%最大误差
meanerrcoif3=mean(errcoif3);%平均误差
stderrcoif3=std(errcoif3);%标准差


plot(P,'.','Color','r','MarkerSize',10);hold on
plot(iwcoif3,'Color',[0 102 153]/255,'LineWidth',1.1)
plot(errcoif3,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend('原始数据','重构曲线','location','northwest','fontsize', 15, 'fontname', '微软雅黑')

% title([num2str(lcoif3(1)),'项coif3小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
% 
% saveas(gcf,[path,'plan',num2str(plan),'coif3',num2str(sum(lxcoif3(1:1+cutcoif3level))),'.fig'])
% print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'coif3',num2str(sum(lxcoif3(1:1+cutcoif3level)))])
% print(gcf,'-dpng',[path,'plan',num2str(plan),'coif3',num2str(sum(lxcoif3(1:1+cutcoif3level)))])
% print(gcf,'-dsvg',[path,'plan',num2str(plan),'coif3',num2str(sum(lxcoif3(1:1+cutcoif3level)))])
%% sym5

figure
tic
% nexttile
levelsym5=3;
cutsym5level=1;
[csym5,lsym5] = wavedec(P(:,1),levelsym5,'sym5');
rsym5=csym5;
rsym5(sum(lsym5(1:1+cutsym5level))+1:end) = 0;
iwsym5 = waverec(rsym5,lsym5,'sym5');

sum(lsym5(1:1+cutsym5level))

elapsedTimesym5=toc;
energysym5=iwsym5.^2;
errsym5=vecnorm(P-Psym5);
maxerrsym5=max(errsym5);%最大误差
meanerrsym5=mean(errsym5);%平均误差
stderrsym5=std(errsym5);%标准差



plot(P,'.','Color','r','MarkerSize',10);hold on
plot(iwsym5,'Color',[0 102 153]/255,'LineWidth',1.1)

set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend('原始数据','重构曲线','fontsize', 15, 'fontname', '微软雅黑')

% title([num2str(lsym5(1)),'项sym5小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
% 
% saveas(gcf,[path,'plan',num2str(plan),'sym5',num2str(sum(lxsym5(1:1+cutsym5level))),'.fig'])
% print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'sym5',num2str(sum(lxsym5(1:1+cutsym5level)))])
% print(gcf,'-dpng',[path,'plan',num2str(plan),'sym5',num2str(sum(lxsym5(1:1+cutsym5level)))])
% print(gcf,'-dsvg',[path,'plan',num2str(plan),'sym5',num2str(sum(lxsym5(1:1+cutsym5level)))])

%% Fourier变换对比试验
tic
figure
nd=256;
z = frdescp1d(nOP);
P_FFT = ifrdescp1d(z,nd);
elapsedTimef=toc;
errf=abs(nOP-P_FFT);
maxerrf=max(errf);%最大误差
meanerrf=mean(errf);%平均误差
stderrf=std(errf);%标准差
MSEf=immse(nOP,P_FFT);%均方误差
CCf=XGD(nOP,P_FFT);
plot(nOP,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(P_FFT,'Color',[0 102 153]/255,'LineWidth',3)
plot(errf,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(nd),'项Fourier变换'],'fontsize', 15, 'fontname', '微软雅黑')
snrf=SNR(nOP,P_FFT);
%% 累积归一化能量图
figure
energyv=cumsum(sort(abs(Lambda1(4:260)).^2,'descend'));
normenergyv=energyv/energyv(end);
energyf=cumsum(sort(abs(z(1:256)).^2,'descend'));
normenergyf=energyf/energyf(end);
energydb4=cumsum(sort(abs(cdb4(1:256)).^2,'descend'));
nomenergydb4=energydb4/energydb4(end);
energydb2=cumsum(sort(abs(cdb2(1:256)).^2,'descend'));
nomenergydb2=energydb2/energydb2(end);
% energystl=[0;cumsum(sort(sortedlam.^2,'descend'))];
% normenergystl=energystl/energystl(end);
hold on
plot(normenergyv(1:256),'LineWidth',3)
plot(normenergyf(1:256),'LineWidth',3)
plot(nomenergydb4(1:256),'LineWidth',3)
plot(nomenergydb2(1:256),'LineWidth',3)
% plot(normenergystl(1:256),'LineWidth',1.5)
legend({'V系统正交变换','Fourier变换','db4小波变换','db2小波变换'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
xlim([0 256])
ylim([0 1.2])
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('累积归一化频谱能量图','fontsize', 15, 'fontname', '微软雅黑')