clear
plan=2;
switch plan
    case 1
        load bird200.mat
    case 2
        load fire500.mat
        ln=5;
    case 3
        load yezi600.mat
        ln=5;
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
        ln=5;
    case 9
        load G-200.mat
        ln=5;
        
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
    case 21
        load hudie2fjy420.mat      
        gpoint=P;
        ln=5;
    case 22
        load star3fjy360.mat
        gpoint=P;
        ln=5;    
end
load(['tlist' num2str(plan) '.mat'])
P=[gpoint;gpoint(1,:)];
% P=gpoint;
% [C,L]=wavedec2(gpoint,4,'DB4');
% C(L(4)+1:end)=0;
% grec=waverec2(C,L,'DB4');
% plot(gpoint(:,1),gpoint(:,2),'.')
% hold on
% plot(grec(:,1),grec(:,2))
[c,l] = wavedec(P(:,1),4,'db4');
c(l(ln)+1:end) = 0;
iwx = waverec(c,l,'db4');
[c,l] = wavedec(P(:,2),4,'db4');
c(l(ln)+1:end) = 0;
iwy = waverec(c,l,'db4');
P_DWT = [iwx,iwy];
iwx=[iwx;iwx(1)];
iwy=[iwy;iwy(1)];
t=linspace(0,1,length(iwx))';
tt=linspace(0,1,4*length(iwx))';
P_DWx=spline(t,iwx,tt);
P_DWy=spline(t,iwy,tt);
figure
if plan<20
load(['point',num2str(plan),'-200','.mat'])
end
plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10)
hold on
plot(P_DWx,P_DWy,'Color',[0 102 153]/255,'LineWidth',1.1)
legend({'原始数据','拟合曲线'},'location','northwest','fontsize', 15, 'fontname', '微软雅黑')
% wucha=vecnorm((P_DWT-P),2,2);
[wuchaW,pjwuchaW]=distanceerror(P,P_DWT);

tzwuchaW=wuchaW(tlist);

pathname='C:\CodeStore\matlab\vfitcurve\data\';
tzwcname=['tzwuchaW',num2str(plan),'.mat'];
save ([pathname,tzwcname],'tzwuchaW','pjwuchaW','wuchaW') 
axis equal
axis off