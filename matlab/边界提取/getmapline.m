close all
clear all
Map=shaperead('C:\CodeStore\matlab\边界提取\hydl.shp');
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
    maplong(i)=length(Map(i).X);
    if maplong(i)>500
    plot(Map(i).X,Map(i).Y) 
    end
end
[smap,mapidx]=sort(maplong,'descend');
figure
idxs=24;
plot(Map(mapidx(idxs)).X,Map(mapidx(idxs)).Y)   

% for i=1:length(P)
%     pause(0.1)
% plot(P(i,1),P(i,2),'.'),hold on;
% end