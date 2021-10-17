clear
plan=5;
switch plan
    case 1
        load bird200.mat
    case 2
        load fire500.mat
    case 3
        load yezi600.mat
    case 4
        load shizi1500.mat
    case 5
        load fenghuang2000.mat
end
load(['tlist' num2str(plan) '.mat'])
P=[gpoint;gpoint(1,:)];
figure
% [C,L]=wavedec2(gpoint,4,'DB4');
% C(L(4)+1:end)=0;
% grec=waverec2(C,L,'DB4');
% plot(gpoint(:,1),gpoint(:,2),'.')
% hold on
% plot(grec(:,1),grec(:,2))
[c,l] = wavedec(P(:,1),4,'db4');
c(l(5)+1:end) = 0;
iwx = waverec(c,l,'db4');
[c,l] = wavedec(P(:,2),4,'db4');
c(l(5)+1:end) = 0;
iwy = waverec(c,l,'db4');
iwx=[iwx;iwx(1)];
iwy=[iwy;iwy(1)];
P_DWT = [iwx,iwy];
plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15)
hold on
plot(iwx,iwy,'Color',[0 102 153]/255,'LineWidth',1.5)
wucha=vecnorm((P_DWT(1:end-2,:)-gpoint),2,2);
pj=mean(wucha);
tzwuchaW=vecnorm((P_DWT(tlist,:)-gpoint(tlist,:)),2,2);
tzwcname=['tzwuchaW',num2str(plan),'.mat'];
save (tzwcname,'tzwuchaW') 
axis equal
axis off