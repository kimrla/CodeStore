%% 画重新参数化的X轴数轴变化示意图

clear
close all
N=2;
x=linspace(0,1,2^(N-1)+1);
y=zeros(size(x));
plot(x,y,'color','k','linewidth', 1.1);


hold on
x1=0.34;
x2=0.6;
plot([x1 x2],0,'.','MarkerSize',20)

ylim([0 0.01])
set(gca,'ytick',[],'xtick',sort([x x1 x2]),'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
