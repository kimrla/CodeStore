clear
plan=19;
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
    case 6
        load hudie3000.mat
     case 7
        load niao21000.mat
    case 8
        load huacao4-1500.mat
    case 9
        load G-200.mat
        
    case 11
        load bird200r.mat
    case 12
        load fire500r.mat
    case 13
        load yezi600r.mat
    case 14
        load shizi1500r.mat
    case 15
        load fenghuang2000r.mat
    case 16
        load hudie3000r.mat
    case 17
        load niao21000r.mat
    case 18
        load huacao4-1500r.mat
    case 19
        load G-200r.mat
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
wucha=vecnorm((P_DWT(1:end-1,:)-P),2,2);
pj=mean(wucha);
tzwuchaW=vecnorm((P_DWT(tlist,:)-P(tlist,:)),2,2);

pathname='C:\CodeStore\matlab\vfitcurve\data\';
tzwcname=['tzwuchaW',num2str(plan),'.mat'];
save ([pathname,tzwcname],'tzwuchaW','pj') 
axis equal
axis off