clear
plan=9;
switch plan
    case 1
        load bird200.mat
        nd=84;
    case 2
        load fire500.mat
        nd=173;
%         nd=256;
    case 3
        load yezi600.mat
        nd=187;
%         nd=256;
    case 4
        load shizi1500.mat
        nd=474;
    case 5
        load fenghuang2000.mat
        nd=741;
    case 6
        load hudie3000.mat
        nd=810;
     case 7
        load niao21000.mat
        nd=310;
    case 8
        load huacao4-1500.mat
        nd=670;
    case 9
        load G-200.mat
%         N=6;
        nd=119;
%         nd=256;
        
    case 11
        load bird200r.mat
        nd=84;
    case 12
        load fire500r.mat
        nd=156;
    case 13
        load yezi600r.mat
        nd=176;
    case 14
        load shizi1500r.mat
        nd=474;
    case 15
        load fenghuang2000.mat
        nd=741;
%         N=3;
    case 16
        load hudie3000r.mat
%         N=10;
        nd=810;
    case 17
        load niao21000r.mat
%         N=7;
        nd=310;
    case 18
        load huacao4-1500r.mat
%         N=10;
        nd=676;
    case 19
        load G-200r.mat
        %         N=6;
        nd=84;
    case 21
        load hudie2fjy420.mat
        nd=189;
        gpoint=P;
    case 22
        load star3fjy360.mat
        nd=105;
        gpoint=P;
end
load(['tlist' num2str(plan) '.mat'])
% P=[gpoint;gpoint(1,:)];
P=gpoint;

z = frdescp(P);
P_FFT = ifrdescp(z,nd);
% z=fft(P(:,1));
% z(nd+1:end,:)=0;
% P_FFT=ifft(z);


P=[P;P(1,:)];
P_FFT=[P_FFT;P_FFT(1,:)];
t=linspace(0,1,length(P_FFT))';
tt=linspace(0,1,4*length(P_FFT))';
P_FFx=spline(t,P_FFT(:,1),tt);
P_FFy=spline(t,P_FFT(:,2),tt);
figure
if plan<20
load(['point',num2str(plan),'-200','.mat'])
end
plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10)
hold on
plot(P_FFx,P_FFy,'Color',[0 102 153]/255,'LineWidth',1.1)
% plot(P_FFT(:,1),P_FFT(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)
legend({'原始数据','拟合曲线'},'location','northwest','fontsize', 15, 'fontname', '微软雅黑')
[wuchaF,pjwuchaF]=distanceerror(P,P_FFT);

tzwuchaF=wuchaF(tlist);
pathname='C:\CodeStore\matlab\vfitcurve\data\';
tzwcname=['tzwuchaF',num2str(plan),'.mat'];
save ([pathname,tzwcname],'tzwuchaF','pjwuchaF','wuchaF')
axis equal
axis off
