function [PP,pp] = WANQU(P)
% clear all
% Map = shaperead('am.shp') ;
% j = 1180;%5825;
% % j = 2;
% P = [Map(j).X;Map(j).Y]';
% P = P(1:end-1,:);
% % P = gd;
% figure
% plot(P(:,1),P(:,2),'--','Color',[0,102,153]/255,'LineWidth',1.5,'DisplayName','原始曲线');hold on
% axis equal
% axis off
lscale = 0.045;
dscale = 0.048;
% tol = 0.005;
% [ps,ix] = dpsimplify(P,tol);
%% 检测拐点
q = 0;
for i = 2:length(P)-1
    if (P(i,1)-P(i-1,1))*(P(i+1,1)-P(i,1))>0&&(P(i,2)-P(i-1,2))*(P(i+1,2)-P(i,2))>0  %% 该点是否为拐点
        q = q+1;
        gd(q,:) = P(i,:);   %% 储存拐点数据
        in(q) = i;      %% 储存拐点位置
    end
end

%% 求弯曲的宏观顶角与微观顶角
for i = 2:length(in)-1  %% 宏观顶角
    AB = [gd(i,1)-gd(i-1,1),gd(i,2)-gd(i-1,2)];
    AC = [gd(i,1)-gd(i+1,1),gd(i,2)-gd(i+1,2)];
    in(2,i) = acosd((AB(1)*AC(1)+AB(2)*AC(2))/((sqrt((AB(1))^2+(AB(2))^2))*sqrt((AC(1))^2+(AC(2))^2)));
end
in(3,:) =1;
for i = 2:length(in)-1
    if in(2,i)>70
        if in(2,i+1)<70
            in(3,i) = 0;
        else
                if in(2,i+1)>in(2,i)
                    in(3,i+1) = 0;
                else 
                    in(3,i) = 0;
                end
        end
    else
        if in(2,i+1)>70
            in(3,i+1) = 0;
        end
    end
end
k = 0;
for i = 1:length(in)
    if in(3,i) == 1
        k = k+1;
        gd1(k,:) = gd(i,:);
        in1(:,k) = in(:,i);
    end
end
for i = 2:length(in1)-1  %% 宏观顶角
    AB = [gd1(i,1)-gd1(i-1,1),gd1(i,2)-gd1(i-1,2)];
    AC = [gd1(i,1)-gd1(i+1,1),gd1(i,2)-gd1(i+1,2)];
    in1(2,i) = acosd((AB(1)*AC(1)+AB(2)*AC(2))/((sqrt((AB(1))^2+(AB(2))^2))*sqrt((AC(1))^2+(AC(2))^2)));
end
% in(3,:) =1;
for i = 2:length(in1)-1
    if in1(2,i)>90
        if in1(2,i+1)<90
            in1(3,i) = 0;
        else
                if in1(2,i+1)>in1(2,i)
                    in1(3,i+1) = 0;
                else 
                    in1(3,i) = 0;
                end
        end
    else
        if in1(2,i+1)>90
            in1(3,i+1) = 0;
        end
    end
end
k = 0;
for i = 1:length(in1)
    if in1(3,i) == 1
        k = k+1;
        gd2(k,:) = gd1(i,:);
        in2(:,k) = in1(:,i);
    end
end
for i = 2:length(in2)-1  %% 宏观顶角
    AB = [gd2(i,1)-gd2(i-1,1),gd2(i,2)-gd2(i-1,2)];
    AC = [gd2(i,1)-gd2(i+1,1),gd2(i,2)-gd2(i+1,2)];
    in2(2,i) = acosd((AB(1)*AC(1)+AB(2)*AC(2))/((sqrt((AB(1))^2+(AB(2))^2))*sqrt((AC(1))^2+(AC(2))^2)));
end
% in(3,:) =1;
for i = 2:length(in2)-1
    if in2(2,i)>90
        if in2(2,i+1)<90
            in2(3,i) = 0;
        else
                if in2(2,i+1)>in1(2,i)
                    in2(3,i+1) = 0;
                else 
                    in2(3,i) = 0;
                end
        end
    else
        if in2(2,i+1)>90
            in2(3,i+1) = 0;
        end
    end
end
k = 0;
for i = 1:length(in2)-2
    if in2(3,i) == 1
        k = k+1;
        gd3(k,:) = gd2(i,:);
        in3(:,k) = in2(:,i);
    end
end
for i = 1:2:length(in2)-2
    A = gd2(i,:);
    B = gd2(i+2,:);
    L = sqrt((B(1)-A(1))^2+(B(2)-A(2))^2);
    if L<lscale
        C(1) = P(in2(1,i+2),2)-P(in2(1,i),2);
        C(2) = P(in2(1,i),1)-P(in2(1,i+2),1);
        C(3) = -C(2)*P(in2(1,i),2)-C(1)*P(in2(1,i),1);
        distance = 0;
        for j = in2(1,i):in2(1,i+2)
            distan = abs((C(1)*P(j,1)+C(2)*P(j,2)+C(3))/(sqrt(C(1)^2+C(2)^2)));
            if distan>distance
                distance = distan;
            end
        end
        if distance<dscale
            P(in2(1,i)+1:in2(1,i+2)-1,:) = 0;
        end
    end
end
h = 0;
for i = 1:length(P)
    if P(i,1) ~= 0
        h =h+1;
        PP(h,:) = P(i,:);
        pp(h) = i;
    end
end
% figure
% plot(ps(:,1),ps(:,2),'.','MarkerFaceColor',[255,102,102]/255,'MarkerSize',20,'DisplayName','特征点');hold on
% % plot(PP(:,1),PP(:,2),'.','MarkerFaceColor',[255,255,102]/255,'MarkerSize',20,'DisplayName','拐点');hold on
% plot(PP(:,1),PP(:,2),'Color',[255,102,102]/255,'LineWidth',1.5,'DisplayName','化简曲线');
% legend
% axis equal
% axis off
% legend
% for i = 2:length(in)-1  %% 微观顶角
%     BA = [P(in(i),1)-P(in(i)-1,1),P(in(i),2)-P(in(i)-1,2)];
%     BC = [P(in(i),1)-P(in(i)+1,1),P(in(i),2)-P(in(i)+1,2)];
%     in(i,2) = acos((BA(1)*BC(1)+BA(2)*BC(2))/(abs(BA(1)*BA(2))*abs(BC(1)*BC(2))));
% end