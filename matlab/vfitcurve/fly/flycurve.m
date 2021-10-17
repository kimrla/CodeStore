clear
plan=5;
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
end
        load(['tlist' num2str(plan) '.mat'])
P=[gpoint;gpoint(1,:)];
z = frdescp(P);
P_FFT = ifrdescp(z,nd);
P_FFT=[P_FFT;P_FFT(1,:)];
plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15)
hold on
plot(P_FFT(:,1),P_FFT(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
wucha=vecnorm((P_FFT(1:end-3,:)-gpoint),2,2);
pj=mean(wucha);
tzwuchaF=vecnorm((P_FFT(tlist,:)-gpoint(tlist,:)),2,2);
tzwcname=['tzwuchaF',num2str(plan),'.mat'];
save (tzwcname,'tzwuchaF') 
axis equal
axis off
