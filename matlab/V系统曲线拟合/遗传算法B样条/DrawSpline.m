function DrawSpline(n, p, Pi, ui,a,b)
% B�����Ļ�ͼ����
% ��֪n+1�����ƶ���P(i), k��B������P��2*(n+1)�������ƶ�������, �ڵ�����NodeVector
% plot(Pi(1:n+1,1), Pi(1:n+1,2),...
%                     'o','LineWidth',1,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','g',...
%                     'MarkerSize',6);
% line(Pi(1:n+1,1), Pi(1:n+1,2));
Njp_u = zeros(1, n+1);
j=0;
for u = a : 0.005 : b-0.005
    j=j+1;
    for i = 1 : n+1
        Njp_u(j, i) = Njp(i, p , u, ui);
    end
    p_u(j,:) = Njp_u(j,:)*Pi;

        tempx(j) = p_u(1,1);
        tempy(j) = p_u(1,2);

end
plot(p_u(:,1),p_u(:,2),'Color',[0 102 153]/255,'LineWidth',3)
end