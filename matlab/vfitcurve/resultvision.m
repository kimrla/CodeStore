% 画不同方法拟合曲线的特征误差对比图
clear
plan=2;
switch plan
%     case 1
%         load bird200.mat
%     case 2
%         load fire500.mat
%     case 3
%         load yezi600.mat
%     case 4
%         load shizi1500.mat
%     case 5
%         load fenghuang2000.mat    
%     case 6
%         load hudie3000.mat
%      case 7
%         load niao21000.mat
%     case 8
%         load huacao4-1500.mat
%     case 9
%         load G-200.mat
    case 1
        load bird200r.mat
    case 2
        load fire500r.mat
    case 3
        load yezi600r.mat
    case 4
        load shizi1500r.mat
    case 5
        load fenghuang2000r.mat
    case 6
        load hudie3000r.mat
     case 7
        load niao21000r.mat
    case 8
        load huacao4-1500r.mat
    case 9
        load G-200r.mat
    case 21
        load hudie2fjy420.mat
        gpoint=P;
    case 22
        load star3fjy360.mat
        gpoint=P;
end
figure
load(['tlist' num2str(plan) '.mat'])
load (['tzwuchaV' num2str(plan) '.mat'])
load (['tzwuchaB' num2str(plan) '.mat'])
load (['tzwuchaP' num2str(plan) '.mat'])
load (['tzwuchaF' num2str(plan) '.mat'])
load (['tzwuchaW' num2str(plan) '.mat'])
p=plot([tzwuchaV tzwuchaB tzwuchaP tzwuchaF tzwuchaW],'markersize',10,'LineWidth',1.1);
% p=plot([tzwuchaV tzwuchaP tzwuchaF tzwuchaW],'markersize',10,'LineWidth',1.1);
p(1).Marker='o';
p(2).Marker='+';
p(3).Marker='*';
p(4).Marker='s';
p(5).Marker='.';
% p(1).Color=[1 0.6 0.8];
% p(2).Color=[0.6 0.8 1];
% p(3).Color=[0.6 1 0.8];
% p(4).Color=[0.8 1 0.6];
% p(5).Color=[1 0.8 0.6];
p(1).Color=[0.96, 0.62, 0.24];
p(2).Color=[0.99, 0.57, 0.59];
p(3).Color=[0.10, 0.67, 0.59];
p(4).Color=[0.68, 0.42, 0.89];
p(5).Color=[0.28, 0.55, 0.86];
% p(1).Color='#994A43';
% p(2).Color='#D9E693';
% p(3).Color='#E6857D';
% p(4).Color='#65A1E6';
% p(5).Color='#4B6F99';
legend('V-系统(本文方法)','基于GA的B样条拟合','PIA','Fourier变换','DB4小波变换')
% legend('V-系统(本文方法)','PIA','Fourier变换','DB4小波变换')
xticks(1:1:length(tlist));
xlabel('特征点序号')
ylabel('误差')
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
box off