function PP = Li_openshaw(P)
% clear all
% Map = shaperead('am.shp') ;
% j = 1180;%5825;
% j = 6304;
% P = [Map(j).X;Map(j).Y]';
% P = P(1:end-1,:);
% tol = 0.0005;
% [ps,ix] = dpsimplify(P,tol);
% figure
% plot(P(:,1),P(:,2),'--','Color',[0,102,153]/255,'LineWidth',1.5,'DisplayName','原始曲线');hold on
% axis equal
% axis off
r = 0.019; %  半径
q = 1;
A=P(1,:);  %  设置初始圆的圆心
L = 1;
distance1 = 0;
M(q,:) = P(1,:);

while q<20||distance1>r %设置终止条件
    q = q+1;
    [M(q,:),L,A] = li_op(A,r,P,L);
     distance1 = real(sqrt((A(1)-P(end,1))^2+(A(2)-P(end,2))^2));
end
M(q+1,:) = P(end,:);
PP = M;
% for i =1:length(P)
%     distance1 = sqrt((A(1)-P(end,1))^2+(A(2)-P(end,2))^2);
% if q>100&&distance1<r
%     break;
% else
%     q = q+1;
%     [M(q,:),i] = li_op(A,r,P,i);
%     A = M(q,:);
% end
% end
% figure
% plot(ps1(:,1),ps1(:,2),'.','MarkerFaceColor',[255,102,102]/255,'MarkerSize',20,'DisplayName','特征点');hold on
% plot(M(:,1),M(:,2),'.','MarkerFaceColor',[255,255,102]/255,'MarkerSize',20,'DisplayName','拐点');hold on
% plot(M(:,1),M(:,2),'Color',[255,102,102]/255,'LineWidth',1.5,'DisplayName','化简曲线');
% axis equal
% axis off
% legend