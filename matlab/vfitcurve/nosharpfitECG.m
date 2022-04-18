clear
clc
close all;
%------ 具体数据------------------------------------------------------
PATH= 'C:\CodeStore\matlab\mit-bih-arrhythmia-database-1.0.0'; % path, 这里就是写刚才你保存的数据地址
ecgname='106'; %101 100+1100  %106 100+1300 %112 140+1100
nOPnum=1300;%截取数据个数
stp=100;%截取数据点起点
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
% title('原始信号','fontsize', 15, 'fontname', '微软雅黑')
hold on
nOP=ORIGINP;

% figure %中值滤波
% R=100;%中值滤波的窗口W=2*R+1
% nOP=myMedfilt(ORIGINP, R);%中值滤波后的结果
% plot(nOP,'LineWidth',3);
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% title('中值滤波处理后的信号','fontsize', 15, 'fontname', '微软雅黑')


% figure
% [c,l] = wavedec(nOP,3,'coif5');%小波分析滤波 3层分解
% a3 = wrcoef('a',c,l,'coif5');
% plot(a3,'LineWidth',3);
% nOP=a3;

% [XC1,CXC1,LXC1,PERF0_1,PERF2_1]=wdenoise(nOP);%小波去噪
% nOP=XC1';
% plot (XC1,'LineWidth',3) %XC1为去噪后的数据
PEAK_1=ATRPOSITION(ATRPOSITION>=stp & ATRPOSITION<stp+nOPnum)-stp+1;
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% title('小波滤波处理后的信号','fontsize', 15, 'fontname', '微软雅黑')
x=linspace(0,1,nOPnum)';
%%
f=nOP;
t=x;
% PEAK_1=getecgtezheng(wdenoise(a3))';
tzloc=[1;PEAK_1;nOPnum];
% plot(tzloc,nOP(tzloc),'o','color','r','MarkerSize',10)

tezhengt=t(tzloc);
k=3;
N=6;

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



path='C:\CodeStore\matlab\vfitcurve\nosharp\';
draw2tree(k,N,Lambda)
drawheat(k,N,Lambda)
%% 本方法最终重构结果
figure
% tiledlayout(1,4,'TileSpacing','compact')
% nexttile
curve=LSMatrix_V(k,N,newt)*Lambda;
err=abs(f-curve);
maxerr=max(err);%最大误差
meanerr=mean(err);%平均误差
stderr=std(err);%标准差
MSE=immse(f,curve);%均方误差
CC=XGD(f,curve);
plot(f,'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
plot(curve,'Color',[0 102 153]/255,'LineWidth',3)
plot(err,'linewidth',3)
yrange=get(gca,'Ylim');
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0],'Ylim',yrange)
% title('最终重构结果','fontsize', 15, 'fontname', '微软雅黑')
% lgd=legend({'原始信号','重构信号','绝对误差'},'fontsize', 15, 'fontname', '微软雅黑');
% lgd.Layout.Tile = 5;
snr1=SNR(f,curve);
saveas(gcf,[path,'ecg',ecgname,'V','.fig'])
print(gcf,'-depsc','-tiff',[path,'ecg',ecgname,'V'])
print(gcf,'-dpng',[path,'ecg',ecgname,'V'])
print(gcf,'-dsvg',[path,'ecg',ecgname,'V'])
% plot(gpoint(locs,1),gpoint(locs,2),'s','MarkerSize',10)
% plot(gpoint(tlist,1),gpoint(tlist,2),'s','MarkerSize',10)
% duandian=LSMatrix_V(k,N,CList(:,1))*Lambda;
% plot(CList(:,1),duandian,'s','MarkerSize',10)
%% 取部分基函数重构结果
% range2cut=1:16;
% [~,slidx]=sort(abs(Lambda),'descend');
% 
% partlambda=Lambda(slidx(range2cut));
% partmat=LSMatrix_V(k,N,newt);
% partmat=partmat(:,slidx(range2cut));
% partcurve=partmat*partlambda;
% figure
% plot(f,'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
% plot(partcurve,'Color',[0 102 153]/255,'LineWidth',3)


%% db4小波变换对比实验
figure
% nexttile
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
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0],'Ylim',yrange)
% title([num2str(ldb4(1)),'项db4小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
snrdb4=SNR(f,iwpdb4);
saveas(gcf,[path,'ecg',ecgname,'db4','.fig'])
print(gcf,'-depsc','-tiff',[path,'ecg',ecgname,'db4'])
print(gcf,'-dpng',[path,'ecg',ecgname,'db4'])
print(gcf,'-dsvg',[path,'ecg',ecgname,'db4'])
%% db2小波变换对比实验
figure
% nexttile
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
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0],'Ylim',yrange)
% title([num2str(ldb2(1)),'项db2小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
snrdb2=SNR(f,iwpdb2);
saveas(gcf,[path,'ecg',ecgname,'db2','.fig'])
print(gcf,'-depsc','-tiff',[path,'ecg',ecgname,'db2'])
print(gcf,'-dpng',[path,'ecg',ecgname,'db2'])
print(gcf,'-dsvg',[path,'ecg',ecgname,'db2'])
%% Fourier变换对比试验
figure
% nexttile
nd=128;
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

set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0],'Ylim',yrange)


% title([num2str(nd),'项Fourier变换'],'fontsize', 15, 'fontname', '微软雅黑')
snrf=SNR(f,P_FFT);
saveas(gcf,[path,'ecg',ecgname,'f','.fig'])
print(gcf,'-depsc','-tiff',[path,'ecg',ecgname,'f'])
print(gcf,'-dpng',[path,'ecg',ecgname,'f'])
print(gcf,'-dsvg',[path,'ecg',ecgname,'f'])

% figure
% plot(f,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
% plot(P_FFT,'Color',[0 102 153]/255,'LineWidth',3)
% plot(errf,'linewidth',3) 
% legend({'原始信号','重构信号','绝对误差'},'location','eastoutside','fontsize', 30, 'fontname', '微软雅黑');
% % legend({'原始信号','重构信号','绝对误差'},'location','eastoutside','fontsize', 10, 'fontname', '微软雅黑');
% saveas(gcf,[path,'ecg',ecgname,'lgd','.fig'])
% print(gcf,'-depsc','-tiff',[path,'ecg',ecgname,'lgd'])
% print(gcf,'-dpng',[path,'ecg',ecgname,'lgd'])
% print(gcf,'-dsvg',[path,'ecg',ecgname,'lgd']) %lgd legend

%% 累积归一化能量图
figure
energyv=cumsum(sort(abs(Lambda(4:128)).^2,'descend'));
normenergyv=energyv/energyv(end);
% energyf=cumsum(sort(abs(z(4:128)).^2,'descend'));
energyf=cumsum(sort(abs(z([length(z)/2-64:length(z)/2-2,length(z)/2+2:length(z)/2+64])).^2,'descend'));
normenergyf=energyf/energyf(end);
energydb4=cumsum(sort(abs(cdb4(4:128)).^2,'descend'));
nomenergydb4=energydb4/energydb4(end);
energydb2=cumsum(sort(abs(cdb2(4:128)).^2,'descend'));
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
xlim([0 128])
ylim([0 1.2])
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0])

saveas(gcf,[path,'ecg',ecgname,'NL','.fig'])
print(gcf,'-depsc','-tiff',[path,'ecg',ecgname,'NL'])
print(gcf,'-dpng',[path,'ecg',ecgname,'NL'])
print(gcf,'-dsvg',[path,'ecg',ecgname,'NL'])
% title('累积归一化频谱能量图','fontsize', 15, 'fontname', '微软雅黑')
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