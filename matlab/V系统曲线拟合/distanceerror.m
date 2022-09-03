function [error,ave,maxda,scha] = distanceerror(P,PP)
%%%该函数为计算化简曲线PP距离P中每个点的最近距离，可将P作为特征点集，PP作为化简曲线。因为化简曲线在实际为折线，因此这个程序
%实际上计算的是特征点距离折线的最短距离，这个距离有可能是到P中点的距离，也有可能是与折线间的距离。
sum = 0; % 初始化距离总误差

for i = 1:length(P)  %确定循环次数，确定距离特征点集P中每个点的最小距离
    q  = 0;  % 记录距离P(i)点最近的折线段端点序号，初始化为0
    distance = 1000; %距离误差初始化为某一最大值，以寻找最小距离
    for j = 2:length(PP)-1  %确定距离P(i)点最近的折线段
        distan = sqrt((PP(j,1)-P(i,1))^2+(PP(j,2)-P(i,2))^2);
        if distance>distan
            distance = distan;
            q = j;  % 按照距离确定最终距离最小的折线段，记录其端点序号
        end
    end
    
 %%
 %求出该点到这个折线的最短距离，因为一个端点连接两个折线，因此要分别求出两个折线的方程并比较两个折线到P（i）的距离，最终结果取距离小的那个
    C(1) = PP(q+1,2)-PP(q,2);  %% 求出该段线方程系数
    C(2) = PP(q,1)-PP(q+1,1);
    C(3) = -C(2)*PP(q,2)-C(1)*PP(q,1);
    distance1 = abs((C(1)*P(i,1)+C(2)*P(i,2)+C(3))/(sqrt(C(1)^2+C(2)^2)));  %% 该折线到圆心的距离
    C1(1) = PP(q,2)-PP(q-1,2);  %% 求出该段线方程系数
    C1(2) = PP(q-1,1)-PP(q,1);
    C1(3) = -C1(2)*PP(q-1,2)-C1(1)*PP(q-1,1);
    distance2 = abs((C1(1)*P(i,1)+C1(2)*P(i,2)+C1(3))/(sqrt(C1(1)^2+C1(2)^2)));  %% 该折线到P(i)的距离
    if distance1<distance2
        error(i) = distance1;
    else
        error(i) = distance2;
    end
    sum = sum+error(i); %累加距离误差
    
end
ave = sum/length(P); %计算平均距离误差
maxda=max(error);%最大误差

sum1=0;
for  i = 1:length(P)
    sum1=sum1+(error(i)- ave)^2;
end
scha=sqrt(   sum1/length(P)    );%标准差

end