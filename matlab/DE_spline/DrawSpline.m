function DrawSpline(n, p, Pi, ui)
% B样条的绘图函数
% 已知n+1个控制顶点P(i), k次B样条，P是2*(n+1)矩阵存控制顶点坐标, 节点向量NodeVector
plot(Pi(1:n+1,1), Pi(1:n+1,2),...
                    'o','LineWidth',1,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',6);
line(Pi(1:n+1,1), Pi(1:n+1,2));
Njp_u = zeros(1, n+1);
j=0;
for u = 0 : 0.005 : 1-0.005
    j=j+1;
    for i = 1 : n+1
        Njp_u(j, i) = Njp(i, p , u, ui);
    end
    p_u = Njp_u(j,:)*Pi;
    if u == 0
        tempx = p_u(1,1);
        tempy = p_u(1,2);
        line([tempx p_u(1,1)], [tempy p_u(1,2)],...
            'Marker','.','LineStyle','-', 'Color',[.3 .6 .9], 'LineWidth',3);
    else
        line([tempx p_u(1,1)], [tempy p_u(1,2)],...
            'Marker','.','LineStyle','-', 'Color',[.3 .6 .9], 'LineWidth',3);
        tempx = p_u(1,1);
        tempy = p_u(1,2);
    end
end