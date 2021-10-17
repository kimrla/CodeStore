%%  下载曲线1,2,3,4 
load P1.mat;
%load P2.mat;
%load P3.mat;
% load P4.mat;
% load P5.mat;
figure
   plot(P(:,1),P(:,2),'r','Linewidth',1);hold on  %展示原海岸线
   axis off

     %%     %%  找海岸线特征点(DP算法)  % 找原始曲线的特征点 ps1
    [ps,ix] = dpsimplify(P,0.0001);   
   %%  Li_openshaw算法
   PP=Li_openshaw(P);
   
   figure
plot(ps(:,1),ps(:,2),'.','MarkerFaceColor',[255,102,102]/255,'MarkerSize',20,'DisplayName','特征点');hold on
plot(PP(:,1),PP(:,2),'.','MarkerFaceColor',[255,255,102]/255,'MarkerSize',20,'DisplayName','拐点');hold on
plot(PP(:,1),PP(:,2),'Color',[255,102,102]/255,'LineWidth',1.5,'DisplayName','化简曲线');
axis equal
axis off