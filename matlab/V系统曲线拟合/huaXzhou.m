%% 画重新参数化的X轴数轴变化示意图

clear
close all
N=4;
x=linspace(0,1,2^(N-1)+1);
% x=[0,0.25,0.375,0.5,1];
x=[0,0.25,0.5,1];
y=zeros(size(x));
plot(x,y,'color','k','linewidth', 1.1);
ylim([0 0.01])
set(gca,'ytick',[],'xtick',x,'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',[0 0 0 0])

% set (gca,'position',[0.1,0.1,0.8,0.1] )
% set (gcf,'unit','normalized','position',[0.1,0.1,0.8,0.1] )
set (gcf,'position',get(gcf,'position').*[1 1 1 0.1])
% box off
% ax1 = axes('Position',get(gca,'Position'),...
%     'XAxisLocation','bottom',...
%     'YAxisLocation','right',...
%    'Color','none',...
%     'XColor','k','YColor','k');
% set(ax1,'XTick',[]);
% set(ax1,'YTick',[]);
% box on;
hold on
x1=0.34;
x2=0.6;
plot([x1 x2],0,'.','MarkerSize',20)

ylim([0 0.01])
set(gca,'ytick',[],'xtick',sort([x x1 x2]),'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
