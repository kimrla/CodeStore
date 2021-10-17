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
Pp=[gpoint;gpoint(1,:)];
k = 3;
n = k + 1;
pro = 1;
NumItr =3;
[PP_all,T] = PIA_CurApr(Pp(:,1:2),k,pro,NumItr);
U = linspace(0,1,10*length(Pp));                  % B样条曲线采样点对应参数值
C1 = bspline_deboor(n,T,[real(PP_all(:,NumItr+1)),imag(PP_all(:,NumItr+1))],U);% 迭代NumItr步后曲线
plot(gpoint(:,1),gpoint(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15)
hold on
plot(C1(:,1),C1(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
Chord = vecnorm(diff(Pp,1,1),2,2); 
normt = [0;cumsum(Chord)/sum(Chord)];  
Cd = bspline_deboor(n,T,[real(PP_all(:,NumItr+1)),imag(PP_all(:,NumItr+1))],normt);
wucha=vecnorm((Cd(1:end-1,:)-gpoint),2,2);
pjwucha=mean(wucha);
tzwuchaP=wucha(tlist);
tzwcname=['tzwuchaP',num2str(plan),'.mat'];
save (tzwcname,'tzwuchaP') 
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