close,clear
Map=shaperead('C:\CodeStore\matlab\边界提取\山脉.shp');
% mapidx=[];
% P=[Map(mapidx).X;Map(mapidx).Y]';
P=[Map.X;Map.Y]';
P(isnan(P(:,1)),:)=[];
P=unique(P,'rows','stable');
% scatter(P(:,1),P(:,2))
% plot(P(:,1),P(:,2))
figure 
hold on
for i=1:length(Map)
    plot(Map(i).X,Map(i).Y)    
end
% for i=1:length(P)
%     pause(0.1)
% plot(P(i,1),P(i,2),'.'),hold on;
% end