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
figure
load(['tlist' num2str(plan) '.mat'])
load (['tzwuchaV' num2str(plan) '.mat'])
load (['tzwuchaB' num2str(plan) '.mat'])
load (['tzwuchaP' num2str(plan) '.mat'])
load (['tzwuchaF' num2str(plan) '.mat'])
load (['tzwuchaW' num2str(plan) '.mat'])
p=plot([tzwuchaV tzwuchaB tzwuchaP tzwuchaF tzwuchaW],'markersize',10);
p(1).Marker='o';
p(2).Marker='+';
p(3).Marker='*';
p(4).Marker='s';
p(5).Marker='.';
legend('V-系统(本文方法)','基于GA的B样条拟合','PIA','Fourier变换','DB4小波变换')
xticks(1:1:length(tlist));
xlabel('特征点序号')
ylabel('误差')
box off