clear
plan=9;
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
legend('V-系统(本文方法)','基于GA的B样条拟合','PIA','Fourier变换','DB4小波变换')
% legend('V-系统(本文方法)','PIA','Fourier变换','DB4小波变换')
xticks(1:1:length(tlist));
xlabel('特征点序号')
ylabel('误差')
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
box off