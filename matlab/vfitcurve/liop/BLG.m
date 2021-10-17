function P1 = BLG(P)
% clear all;
% Map = shaperead('am.shp') ;
% j = 3088;%5825;
% % j = 2;
% P = [Map(j).X;Map(j).Y]';
% P = P(1:end-1,:);
max = 0;
if P(1,:) ==P(end,:)
    for i = 1:length(P)
        distance = sqrt((P(1,1)-P(i,1))^2+(P(1,2)-P(i,2))^2);
        if distance>max
            max = distance;
            ip = i;
        end
    end
    P1 = P;
    up = 10;
    P1(1,4) = 10;
    P1(ip,4) = 10;
    P1(end,4) = 10;
    q = 0;
    P1 = BLGTree(1,ip,P1,q,up);
    P1 = BLGTree(ip,length(P),P1,q,up);
else
    P1 = P;
    up =10;
    P1(1,4) = 10;
    P1(end,4) = 10;
    q = 0;
    P1 = BLGTree(1,length(P),P1,q,up);
end
for i = 1:length(P1)
    P1(i,5) = i;
end
a4 = P1(:,4);
[a4,pos]=sort(a4,"descend");
a1 = P1(pos,1);
a2=P1(pos,2);
a3=P1(pos,3);%a2,a3是排列好的第二行和第三行
a5 = P1(pos,5);
P1 = [a1,a2,a3,a4,a5];

% [a3,pos]=sort(a3);
% a1 = P1(pos,1);
% a2=P1(pos,2);
% a4=P1(pos,4);%a2,a3是排列好的第二行和第三行
% a5 = P1(pos,5);
% P1 = [a1,a2,a3,a4,a5];