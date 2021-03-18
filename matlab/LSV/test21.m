
%% 例1:对圆采样，用1次V做带约束最小二乘逼近

clear all
close all

%% Test1
t = [1/64:1/32:1]';
x = cos(2*pi*t);
y = sin(2*pi*t);
% figure,
% plot(x,y,'.','MarkerSize',10);
% axis equal

k = 1;
N = 4;

%% 约束方程
[DR,DL] = VContinuityInfo1(N);
% CList = [0,1,0;
%          0.5,0.5,0;
%          0.25,0.25,0;
%          0.75,0.75,0];
CList = [0.5,0.5,0;
         0.25,0.25,0;
         0.75,0.75,0;
         1/8,1/8,0;
         3/8,3/8,0;
         5/8,5/8,0;
         7/8,7/8,0;
         0,1,0];
p = sum(CList(:,3)+1);  % 约束方程数目
SegNum = 2^(N-1); % 分段数
VNum = 2^N;       % 基函数数目
knots = linspace(0,1,SegNum+1);% 节点向量
C = zeros(p,VNum);    % 约束矩阵
csidx = 0; 
for h = 1 : length(CList)           % 给每个约束建立相应的方程组
    x1 = find(knots==CList(h,1));   % 左节点
    x2 = find(knots==CList(h,2));   % 右节点
    for c = 0 : CList(h,3)          % 从C^0约束至C^CList(h,3)约束，每个约束建立一个方程
        csidx = csidx + 1;
        C(csidx,:) = [DR(:,x1,c+1) - DL(:,x2,c+1)]';
    end
end  

%% 系数矩阵
A = LSMatrix_V(k,N,t);

M = [2*A'*A,C';
     C,zeros(p)];
d = zeros(p,1);
bx = [2*A'*x;d]; 
by = [2*A'*y;d]; 

X = M\bx;
Y = M\by;

X = X(1:VNum);
Y = Y(1:VNum);

%% Show
InfoV1Bas = BaseGene_V1(N);   % 1次V系统基函数信息(非离散采样)
[VRInfo,NumSeg] = VReconstruction_Polyline([X Y],InfoV1Bas);
figure,
plot(x,y,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
for j = 1 : NumSeg
    plot(VRInfo(j,5:6),VRInfo(j,7:8),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
end
axis equal
axis off