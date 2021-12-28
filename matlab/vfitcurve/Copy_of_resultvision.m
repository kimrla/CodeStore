clear
plist=[2,3,8,9,21,22];
% plist=[9];
for i=1:length(plist)

plan=plist(i);
load(['tlist' num2str(plan) '.mat'])
load (['tzwuchaV' num2str(plan) '.mat'])
load (['tzwuchaB' num2str(plan) '.mat'])
% load (['tzwuchaP' num2str(plan) '.mat'])
load (['tzwuchaF' num2str(plan) '.mat'])
load (['tzwuchaW' num2str(plan) '.mat'])
% pjtzwcV(i,1)=mean(tzwuchaV);
% pjtzwcB(i,1)=mean(tzwuchaB);
% pjtzwcP(i,1)=mean(tzwuchaP);
% pjtzwcF(i,1)=mean(tzwuchaF);
% pjtzwcW(i,1)=mean(tzwuchaW);
ymax=max([wuchaV wuchaB wuchaF wuchaW]);
figure
% subplot(4,1,1),plot(wuchaV','markersize',10,'LineWidth',1.1);
% subplot(4,1,2),plot(wuchaB','markersize',10,'LineWidth',1.1);
% subplot(4,1,3),plot(wuchaF','markersize',10,'LineWidth',1.1);
% subplot(4,1,4),plot(wuchaW','markersize',10,'LineWidth',1.1);
% t=tiledlayout(4,1);
% nexttile
% plot(wuchaV','markersize',10,'LineWidth',1.1);
% xlim([0 length(wuchaV)]);
% ylim([0 ymax])
% title('V-系统(本文方法)')
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% nexttile
% plot(wuchaB','markersize',10,'LineWidth',1.1);
% xlim([0 length(wuchaV)]);
% ylim([0 ymax])
% title('基于遗传算法的B样条曲线拟合')
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% nexttile
% plot(wuchaF','markersize',10,'LineWidth',1.1);
% xlim([0 length(wuchaV)]);
% ylim([0 ymax])
% title('Fourier变换')
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% nexttile
% plot(wuchaW','markersize',10,'LineWidth',1.1);
% xlim([0 length(wuchaV)]);
% ylim([0 ymax])
% title('DB4小波变换')
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% ylabel(t,'误差','fontsize', 10, 'fontname', '微软雅黑')
% set(t, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% for i=1:4
%     plot3(i,)
% end

%3维图
x=linspace(1,length(wuchaV),length(wuchaV));
y=ones(size(wuchaV));
% [X,Y]=meshgrid(x,[1 2 3 4]);
% z=[wuchaV;wuchaB;wuchaF;wuchaW];
% mesh(X,Y,z);
hold on
box on
% grid on
plot3(x,1*y,wuchaV,'LineWidth',1.1);
% text(x(tlist),1*y(tlist),wuchaV(tlist),'o','fontsize',10)
plot3(x,2*y,wuchaB,'LineWidth',1.1);
% text(x(tlist),2*y(tlist),wuchaB(tlist),'o','fontsize',10)
plot3(x,3*y,wuchaF,'LineWidth',1.1);
% text(x(tlist),3*y(tlist),wuchaF(tlist),'o','fontsize',10)
plot3(x,4*y,wuchaW,'LineWidth',1.1);
% text(x(tlist),4*y(tlist),wuchaW(tlist),'o','fontsize',10)
for i=1:length(tlist)
    plot3(x(tlist(i))*y,linspace(0,5,length(wuchaV)),0*y,'k','linewidth',1.5);
end
% plot3(x(tlist),4*y(tlist),0*y(tlist));
view(11,17)
xlim([1 length(wuchaV)]);
ylim([0.8 4.2]);
set(gca,'ytick',[1: 4],'yticklabel',{'V系统','遗传算法','Fourier','DB4小波变换'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
zlabel('误差')


% set(gcf,'ylim',[0 ymax]) 
% set(gcf, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% p=plot(,'markersize',10,'LineWidth',1.1);
% legend('V-系统(本文方法)','基于GA的B样条拟合','Fourier变换','DB4小波变换')
% legend('V-系统(本文方法)','PIA','Fourier变换','DB4小波变换')
% xlim([0 length(wuchaV)]);
% % ylim([-0.2 1]);
% % xticks(1:100:length(wuchaV));
% % xlabel('特征点序号')
% ylabel('误差')
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% box off
end
% p=plot([pjtzwcV pjtzwcB pjtzwcP pjtzwcF pjtzwcW],'markersize',10,'LineWidth',1.1);
% p=plot([tzwuchaV tzwuchaP tzwuchaF tzwuchaW],'markersize',10,'LineWidth',1.1);
% p=plot([wuchaV' wuchaB' wuchaF' wuchaW'],'markersize',10,'LineWidth',1.1);
% p(1).Marker='o';
% p(2).Marker='+';
% p(3).Marker='*';
% p(4).Marker='s';
% p(5).Marker='.';
% p(1).LineStyle='-';
% p(2).LineStyle='--';
% p(3).LineStyle=':';
% p(4).LineStyle='-.';
% p(1).Color=[1 0.6 0.8];
% p(2).Color=[0.6 0.8 1];
% p(3).Color=[0.6 1 0.8];
% p(4).Color=[0.8 1 0.6];
% p(5).Color=[1 0.8 0.6];
% p(1).Color=[0.96, 0.62, 0.24];
% p(2).Color=[0.99, 0.57, 0.59];
% p(3).Color=[0.10, 0.67, 0.59];
% p(4).Color=[0.68, 0.42, 0.89];
% p(5).Color=[0.28, 0.55, 0.86];
% p(1).Color='#994A43';
% p(2).Color='#D9E693';
% p(3).Color='#E6857D';
% p(4).Color='#65A1E6';
% p(5).Color='#4B6F99';
% legend('V-系统(本文方法)','基于GA的B样条拟合','Fourier变换','DB4小波变换')
% % legend('V-系统(本文方法)','PIA','Fourier变换','DB4小波变换')
% xlim([0 length(wuchaV)]);
% ylim([-0.2 1]);
% % xticks(1:100:length(wuchaV));
% % xlabel('特征点序号')
% ylabel('误差')
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% box off