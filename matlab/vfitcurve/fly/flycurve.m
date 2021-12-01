clear
plan=19;
switch plan
    case 1
        load bird200.mat
        nd=84;
    case 2
        load fire500.mat
        nd=156;
    case 3
        load yezi600.mat
        nd=176;
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
        nd=676;
    case 9
        load G-200.mat
%         N=6;
        nd=84;
        
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

end
load(['tlist' num2str(plan) '.mat'])
P=[gpoint;gpoint(1,:)];
z = frdescp(P);
P_FFT = ifrdescp(z,nd);
P_FFT=[P_FFT;P_FFT(1,:)];
plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15)
hold on
plot(P_FFT(:,1),P_FFT(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
wucha=vecnorm((P_FFT(1:end-2,:)-P),2,2);
pj=mean(wucha);
tzwuchaF=vecnorm((P_FFT(tlist,:)-P(tlist,:)),2,2);

pathname='C:\CodeStore\matlab\vfitcurve\data\';
tzwcname=['tzwuchaF',num2str(plan),'.mat'];
save ([pathname,tzwcname],'tzwuchaF','pj')
axis equal
axis off
