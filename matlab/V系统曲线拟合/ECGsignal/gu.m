close all
clear
load GF3MtxOrth.mat % 67*129 131*256 131*512 131*1024 131*1025
ecgname='210_512p';
load (['No',ecgname,'.mat'])
plot(PP,'.','Color',[255 102 102]/255,'MarkerSize',15); hold on
GF3mtx=GF3MtxOrth{3};
F=GF3mtx(:,1:end-1);
lambda=F*PP;
curve=F'*lambda;
err=abs(PP-curve);
maxerr=max(err);%最大误差
meanerr=mean(err);%平均误差
stderr=std(err);%标准差
MSE=immse(PP,curve);%均方误差
plot(curve,'Color',[0 102 153]/255,'LineWidth',3)
plot(err,'linewidth',3)
yrange=get(gca,'Ylim');
% legend({'原始信号','重构信号','绝对误差'},'location','best','fontsize', 10, 'fontname', '微软雅黑');
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0],'Ylim',yrange)
path='C:\CodeStore\matlab\V系统曲线拟合\ECGsignal\';
% saveas(gcf,[path,'ecg',ecgname,'GF','.fig'])
% print(gcf,'-depsc','-tiff',[path,'ecg',ecgname,'GF'])
print(gcf,'-dpng',[path,'ecg',ecgname,'GF'])
% print(gcf,'-dsvg',[path,'ecg',ecgname,'GF'])
%% db4小波变换对比实验
tic
figure
% nexttile
[cdb4,ldb4] = wavedec(PP,3,'db4');
rdb4=cdb4;
rdb4(134:end) = 0;
iwpdb4 = waverec(rdb4,ldb4,'db4');
errdb4=abs(PP-iwpdb4);
maxerrdb4=max(errdb4);%最大误差
meanerrdb4=mean(errdb4);%平均误差
stderrdb4=std(errdb4);%标准差
MSEdb4=immse(PP,iwpdb4);%均方误差
% CCdb4=XGD(PP,iwpdb4);
elapsedTimedb4=toc;
plot(PP,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(iwpdb4,'Color',[0 102 153]/255,'LineWidth',3)
plot(errdb4,'linewidth',3)
% legend({'原始信号','重构信号','绝对误差'},'location','best','fontsize', 10, 'fontname', '微软雅黑');
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0],'Ylim',yrange)
% title([num2str(ldb4(1)),'项db4小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
% snrdb4=SNR(PP,iwpdb4);
% saveas(gcf,[path,'ecg',ecgname,'db4','.fig'])
% print(gcf,'-depsc','-tiff',[path,'ecg',ecgname,'db4'])
print(gcf,'-dpng',[path,'ecg',ecgname,'db4'])
% print(gcf,'-dsvg',[path,'ecg',ecgname,'db4'])
%% Fourier变换对比试验
tic
figure
% nexttile
nd=131;
z = frdescp1d(PP);
P_FFT = ifrdescp1d(z,nd);
errf=abs(PP-P_FFT);
maxerrf=max(errf);%最大误差
meanerrf=mean(errf);%平均误差
stderrf=std(errf);%标准差
MSEf=immse(PP,P_FFT);%均方误差
% CCf=XGD(PP,P_FFT);
elapsedTimef=toc;
plot(PP,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(P_FFT,'Color',[0 102 153]/255,'LineWidth',3)
plot(errf,'linewidth',3) 
% legend({'原始信号','重构信号','绝对误差'},'location','best','fontsize', 10, 'fontname', '微软雅黑');
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0],'Ylim',yrange)


% title([num2str(nd),'项Fourier变换'],'fontsize', 15, 'fontname', '微软雅黑')
% snrf=SNR(PP,P_FFT);
% saveas(gcf,[path,'ecg',ecgname,'Fourier','.fig'])
% print(gcf,'-depsc','-tiff',[path,'ecg',ecgname,'Fourier'])
print(gcf,'-dpng',[path,'ecg',ecgname,'Fourier'])
% print(gcf,'-dsvg',[path,'ecg',ecgname,'Fourier'])

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
energyv=cumsum(sort(abs(lambda(4:128)).^2,'descend'));
normenergyv=energyv/energyv(end);
% energyf=cumsum(sort(abs(z(4:128)).^2,'descend'));
energyf=cumsum(sort(abs(z([length(z)/2-64:length(z)/2-2,length(z)/2+2:length(z)/2+64])).^2,'descend'));
normenergyf=energyf/energyf(end);
energydb4=cumsum(sort(abs(cdb4(4:128)).^2,'descend'));
nomenergydb4=energydb4/energydb4(end);
% energydb2=cumsum(sort(abs(cdb2(4:128)).^2,'descend'));
% nomenergydb2=energydb2/energydb2(end);
% energystl=[0;cumsum(sort(sortedlam.^2,'descend'))];
% normenergystl=energystl/energystl(end);
hold on
plot(normenergyv,'LineWidth',3)
plot(normenergyf,'LineWidth',3)
plot(nomenergydb4,'LineWidth',3)
% plot(nomenergydb2,'LineWidth',3)
% plot(normenergystl(1:256),'LineWidth',1.5)
legend({'GF系统','Fourier变换','db4小波变换'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
xlim([0 128])
ylim([0 1.2])
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0])

% saveas(gcf,[path,'ecg',ecgname,'NL','.fig'])
% print(gcf,'-depsc','-tiff',[path,'ecg',ecgname,'NL'])
print(gcf,'-dpng',[path,'ecg',ecgname,'NL'])
save (['data',ecgname,'.mat'])
% print(gcf,'-dsvg',[path,'ecg',ecgname,'NL'])
% title('累积归一化频谱能量图','fontsize', 15, 'fontname', '微软雅黑')