
%% 例1:对圆采样，用1次V做最小二乘逼近

clear all
close all

%% Test1
t = [1/64:1/32:1]';
x = cos(2*pi*t);
y = sin(2*pi*t);
figure,
plot(x,y,'.','MarkerSize',10);
axis equal

k = 1;
N = 4;
A = LSMatrix_V(k,N,t);

X = A\x;
Y = A\y;

%% Show1
InfoV1Bas = BaseGene_V1(N);   % 1次V系统基函数信息(非离散采样)
[VRInfo,NumSeg] = VReconstruction_Polyline([X Y],InfoV1Bas);
figure,
plot(x,y,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
for j = 1 : NumSeg
    plot(VRInfo(j,5:6),VRInfo(j,7:8),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
end
axis equal
axis off
