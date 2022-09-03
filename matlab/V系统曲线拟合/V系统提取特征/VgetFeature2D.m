%% 用V系统特征提取2D平面曲线
close all
clear all
example=1;
switch example
    case 1        
        load fire500.mat        
end
P=gpoint;
k=3;
N=floor(log2(length(P)/(k+1)))+1;
t=linspace(0,1,length(P))';
load tlist2.mat
Lambda=LSMatrix_V(k,N,t)\P;
LambdaA=Lambda;
LambdaA(1:length(Lambda)/2,:)=0;
xijie=LSMatrix_V(k,N,t)*LambdaA;
figure
plot(P(:,1),P(:,2))
figure
plot(xijie(:,1))
figure
plot(xijie(:,2))

figure
plot(xijie(:,1),xijie(:,2),'.')



xijiehe=xijie(:,1)+xijie(:,2);
figure
plot(abs(xijiehe))