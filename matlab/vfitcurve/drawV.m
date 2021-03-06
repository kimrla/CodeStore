function drawV(N, k, lamda)
% B样条的绘图函数
% 已知n+1个控制顶点P(i), k次B样条，P是2*(n+1)矩阵存控制顶点坐标, 节点向量NodeVector
% plot(Pi(1:n+1,1), Pi(1:n+1,2),...
%                     'o','LineWidth',1,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','g',...
%                     'MarkerSize',6);
% line(Pi(1:n+1,1), Pi(1:n+1,2));
V = zeros(1, N);
m=0;
for x = 0: 0.001 : 1-0.001
    m=m+1;
    for n = 1 : N
        if n<=k+1
            group=1;
            i=n;
            j=1;
        else
        group=ceil(log2(n/(k+1))+1);
        i=ceil((n-(k+1)*2^(group-2))/2^(group-2));
        j=mod((n-(k+1)*2^(group-2)),2^(group-2));
        if j==0
            j=2^(group-2);
        end
        end
        V(m, n) = vknij(k,group,i,j,x);
    end
    v_x = V(m,:)*lamda;
    if x == 0
        tempx = v_x(1,1);
        tempy = v_x(1,2);
        line([tempx v_x(1,1)], [tempy v_x(1,2)],...
            'Marker','.','LineStyle','-', 'Color',[.3 .6 .9], 'LineWidth',3);
    else
        line([tempx v_x(1,1)], [tempy v_x(1,2)],...
            'Marker','.','LineStyle','-', 'Color',[.3 .6 .9], 'LineWidth',3);
        tempx = v_x(1,1);
        tempy = v_x(1,2);
    end
end