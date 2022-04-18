clear
plan=3;
switch plan
    case 1
        load bird200.mat
        fenbianlv=3;
    case 2
        load fire500.mat
        fenbianlv=0.35;
    case 3
        load yezi600.mat
        fenbianlv=0.3;
    case 4
        load shizi1500.mat
        fenbianlv=0.1;
    case 5
        load fenghuang2000.mat
        fenbianlv=0.1;
    case 6
        load hudie3000.mat
        fenbianlv=0.2;
     case 7
        load niao21000.mat
        fenbianlv=0.04;
    case 8
        load huacao4-1500.mat
        fenbianlv=0.04;
    case 9
        load G-200.mat
        fenbianlv=0.25;
    case 11
        load bird200r.mat
        fenbianlv=3;
    case 12
        load fire500r.mat
        fenbianlv=2;
    case 13
        load yezi600r.mat
        fenbianlv=1;
    case 14
        load shizi1500r.mat
        fenbianlv=0.1;
    case 15
        load fenghuang2000r.mat
        fenbianlv=0.1;
    case 16
        load hudie3000r.mat
        fenbianlv=0.2;
    case 17
        load niao21000r.mat
        fenbianlv=0.5;
    case 18
        load huacao4-1500r.mat
        fenbianlv=0.05;
    case 19
        load G-200r.mat
        fenbianlv=1;
    case 21
        load hudie2fjy420.mat
        fenbianlv=0.5;
        gpoint=P;
    case 22
        load star3fjy360.mat
        gpoint=P;
        fenbianlv=0.032;
    case 23
        load star3fjy72.mat
        gpoint=P;
        fenbianlv=0.0001;
end
load(['tlist' num2str(plan) '.mat'])
P=[gpoint;gpoint(1,:)];
[ps,ix] = dpsimplify(P,fenbianlv);   %%0.0001
if plan<20
load(['point',num2str(plan),'-200','.mat'])
end
plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10)
% plot(gpoint(1:3:end,1),gpoint(1:3:end,2),'.','Color','r','MarkerSize',10)
hold on
% plot(ps(:,1),ps(:,2),'.','Color',[224 222 58]/255,'MarkerSize',15)

Pp=[ps;gpoint(1,:)];
% P=[gpoint;gpoint(1,:)];
% Pp=P(tlist,:);

% Pp=ps;
k = 3;
n = k + 1;
pro = 1;
NumItr =3;
[PP_all,T] = PIA_CurApr(Pp(:,1:2),k,pro,NumItr);
U = linspace(0,1,10*length(gpoint));                  % B样条曲线采样点对应参数值
C1 = bspline_deboor(n,T,[real(PP_all(:,NumItr+1)),imag(PP_all(:,NumItr+1))],U);% 迭代NumItr步后曲线

plot(C1(:,1),C1(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)
legend({'原始数据','拟合曲线'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
Chord = vecnorm(diff([gpoint;gpoint(1,:)],1,1),2,2);
normt = [0;cumsum(Chord)/sum(Chord)];
Cd = bspline_deboor(n,T,[real(PP_all(:,NumItr+1)),imag(PP_all(:,NumItr+1))],normt);


[wuchaP,pjwuchaP,maxwuchaP,schaP]=distanceerror(P,Cd);
tzwuchaP=wuchaP(tlist);
pathname='C:\CodeStore\matlab\vfitcurve\data\';
tzwcname=['tzwuchaP',num2str(plan),'.mat'];
save ([pathname,tzwcname],'tzwuchaP','pjwuchaP','wuchaP','maxwuchaP','schaP')
% plot(Cd(:,1),Cd(:,2),'.','Color','g','LineWidth',1.5)
% tzwcp=vecnorm((C1(tlist,:)-gpoint(tlist,:)),2,2);
% for i=1:length(gpoint)
% %     %     wucha(i)=min(norm(gpoint(1,:)-C1));
% %     wc(:,1)=gpoint(i,1)-C1(:,1);
% %     wc(:,2)=gpoint(i,2)-C1(:,2);
%     wucha(i)=min(sqrt((C1(:,1)-gpoint(i,1)).^2+(C1(:,2)-gpoint(i,2)).^2));
% end
% pj=mean(wucha);
axis equal
axis off