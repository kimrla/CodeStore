% 重新画遗传算法B样条的拟合结果示意图
for gen=1:GM
    if gen<=length(tempui)        
        figure
        plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10)
        hold on
        DrawSplineB(length(tempui{gen})-p-2,p,tempP{gen},tempui{gen},a,b);
        hold off
        xlim([0 1])
        ylim([-110 110])
        box off
%         legend({'原始数据','拟合曲线'},'location','southeast','fontsize', 15, 'fontname', '微软雅黑')
        set(gca, 'linewidth', 1.3, 'fontsize', 15, 'fontname', '微软雅黑')
        axes('position',[0.6 0.2 0.3 0.3])
        u=0:0.005:1-0.005;
        for j=1:length(tempP{gen})
            for i=1:length(u)
                Njp_u(i,j) = Njp(j, p , u(i), tempui{gen});
            end
            %                 subplot(5,5,j),plot(u,Njp_u(:,j));
            plot(u,Njp_u(:,j),'k','LineWidth',2),hold on;
        end
        clear Njp_u
        set(gca, 'linewidth', 1.3, 'fontsize', 15, 'fontname', '微软雅黑')
        
        
        
        
        fig1name=['curve12Bhb',num2str(gen)];
        savefig(fig1name)
        saveas(gcf,fig1name,'png')
        saveas(gcf,fig1name,'emf')
%         
%         figure
%         u=0:0.005:1-0.005;
%         for j=1:length(tempP{gen})
%             for i=1:length(u)
%                 Njp_u(i,j) = Njp(j, p , u(i), tempui{gen});
%             end
%             %                 subplot(5,5,j),plot(u,Njp_u(:,j));
%             plot(u,Njp_u(:,j),'k','LineWidth',1.3),hold on;
%         end
%         clear Njp_u
%         
%         set(gca, 'linewidth', 1.3, 'fontsize', 15, 'fontname', '微软雅黑')
%         fig2name=['curve12B',num2str(gen),'b'];
%         savefig(fig2name)
%         saveas(gcf,fig2name,'png')
%         saveas(gcf,fig2name,'emf')
    end
end
figure
plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10)
hold on
DrawSplineB(neijiedianshuliang(end)+p,p,bestP,bestui,a,b);
xlim([0 1])
ylim([-110 110])
box off
% legend({'原始数据','拟合曲线'},'location','southeast','fontsize', 15, 'fontname', '微软雅黑')
set(gca, 'linewidth', 1.3, 'fontsize', 15, 'fontname', '微软雅黑')
axes('position',[0.6 0.2 0.3 0.3])
u=0:0.005:1-0.005;
for j=1:length(bestP)
    for i=1:length(u)
        Njp_u(i,j) = Njp(j, p , u(i), bestui);
    end
    %                 subplot(5,5,j),plot(u,Njp_u(:,j));
    plot(u,Njp_u(:,j),'k','LineWidth',2),hold on;
end
clear Njp_u
set(gca, 'linewidth', 1.3, 'fontsize', 15, 'fontname', '微软雅黑')



fig1name=['curve12Bhb',num2str(gen)];
savefig(fig1name)
saveas(gcf,fig1name,'png')
saveas(gcf,fig1name,'emf')
% 
% figure
% u=0:0.005:1-0.005;
% for j=1:length(bestP)
%     for i=1:length(u)
%         Njp_u(i,j) = Njp(j, p , u(i), bestui);
%     end
%     %                 subplot(5,5,j),plot(u,Njp_u(:,j));
%     plot(u,Njp_u(:,j),'k','LineWidth',1.3),hold on;
% end
% 
% set(gca, 'linewidth', 1.3, 'fontsize', 15, 'fontname', '微软雅黑')
% fig2name=['curve12B',num2str(gen),'b'];
% savefig(fig2name)
% saveas(gcf,fig2name,'png')
% saveas(gcf,fig2name,'emf')