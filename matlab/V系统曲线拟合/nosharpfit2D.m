%提取特征但不区分尖锐，遍历连续性的重构方法
close all
clear
plan=21;
switch plan    
    case 1
        load bird200.mat
        N=4;
        rd=3;
    case 2
        load fire500.mat
        N=4;
        rd=2;
        fenbianlv=0.35;
        gpoint(end+1,:)=gpoint(1,:);
    case 3
        load yezi600.mat
        N=4;
        rd=1;
        gpoint(end+1,:)=gpoint(1,:);
        fenbianlv=0.3;
    case 4
        load shizi1500.mat
        N=4;
        rd=0.1;
    case 5
        load fenghuang2000.mat
        N=10;
%         N=3;
        rd=0.01;
    case 6
        load hudie3000.mat
%         N=10;
        N=3;
        rd=0.01;
    case 7
        load niao21000.mat
%         N=7;
        N=3;
        rd=0.2;
        fenbianlv=0.04;
    case 8
        load huacao4-1500.mat
        
%         N=10;
        N=3;
        rd=0.1;
        fenbianlv=0.04;
    case 9
        load G-200.mat
%         N=6;
        N=7;
        rd=1;
        fenbianlv=0.15;
        gpoint(end+1,:)=gpoint(1,:);
    case 21
        load hudie2fjy420.mat
        
        gpoint=P;
        gpoint(end+1,:)=gpoint(1,:);
        N=2;
        fenbianlv=0.5;
    case 22
        load yezishu-1024.mat
        gpoint=P{1,1};
        N=2;
        fenbianlv=0.06;
    case 23
        load baolong-1024.mat
        gpoint=P{1,1};
        fenbianlv=0.2;
        N=2;
    case 24
        load sanjiaolongfjy3260.mat
        gpoint=P;
        fenbianlv=0.02;
        N=2;
    case 25
        load bat-1024.mat
        gpoint=P{1,1};
        fenbianlv=0.1;
        N=2;
    case 26
        load jianlong-1024.mat
        gpoint=P{1,1};
        fenbianlv=0.02;
        N=2;
end
runtime=30;
tic
P=gpoint;
num=length(P);
k=3;
t=canshuhua(P);
% t=linspace(0,1,num)';
for i=2:num-1
    [phi(i),K(i)]=qulv(gpoint,i,num);
%     [K(i),~] = PJcurvature(gpoint,i);
    
end
K(num)=K(1);
q=zeros(num,1);
for i=1:num
    r=i+1;
    l=i-1;
    if i==1
        l=num-1;
    elseif i==num
        r=2;
    end
%     mK(i)=1.5*mean([K(l),K(i),K(r)]);
    if K(i)>K(r) && K(i)>K(l)  && K(i)>mean(K)
        %         && phi(i)>pi/6        &&K(i)>3*mean(K)    K(i)>mean(K) phi(i)>pi/6 
        q(i)=1;
    end
    
end

q(1)=1;
q(num)=1;
tzloc=find(q==1);

figure
hold on
plot(gpoint(tzloc,1),gpoint(tzloc,2),'gs','color',[0 102 153]/255,'MarkerSize',10);
% load(['point',num2str(plan),'-200','.mat']);
% load yezi200.mat
plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10,'linewidth',1.5);
legend('原始数据','特征点','location','northeast','fontsize', 15, 'fontname', '微软雅黑')
axis equal
axis off

tezhengt=t(tzloc);
% tezhengt=linspace(0,1,length(tzloc))';
for i=1:length(tezhengt)
    newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),N-1));
end
while length(newtezhengt)-length(unique(newtezhengt))
    N=N+1;
    for i=1:length(tezhengt)
        newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),N-1));
    end
end

newt(tzloc)=newtezhengt;
for i=2:length(tzloc)
    j=tzloc(i-1)+1:tzloc(i)-1;
    newt(j)=(t(j)-tezhengt(i-1))*(newtezhengt(i)-newtezhengt(i-1))/(tezhengt(i)-tezhengt(i-1))+newtezhengt(i-1);
end

newt=newt';
% newt=t;




% if (k+1)*2^(N-1)>length(P)
%     tt=unique([linspace(0,1,(k+1)*2^(N-1))';t]);
%     P=interp1(t,P,tt,'linear');
%     t=tt;
%     num=length(P);
%     for i=2:num-1
%         phi(i)=qulv(P,i,num);
%         [K(i),~] = PJcurvature(P,i);
%     end
%     K(num)=K(1);
%     q=zeros(num,1);
%     for i=1:num
%         r=i+1;
%         l=i-1;
%         if i==1
%             l=num-1;
%         elseif i==num
%             r=2;
%         end
%         if K(i)>K(r) && K(i)>K(l)  && phi(i)>pi/6
%             %         && phi(i)>pi/6        &&K(i)>3*mean(K)
%             q(i)=1;
%         end
%         
%     end
%     q(1)=1;
%     q(num)=1;
%     tzloc=find(q==1);
%     tezhengt=t(tzloc);
%     
%     for i=1:length(tezhengt)
%         newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),N-1));
%     end
%     while length(newtezhengt)-length(unique(newtezhengt))
%         N=N+1;
%         for i=1:length(tezhengt)
%             newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),N-1));
%         end
%     end
%     
%     newt(tzloc)=newtezhengt;
%     for i=2:length(tzloc)
%         j=tzloc(i-1)+1:tzloc(i)-1;
%         newt(j)=(t(j)-t(tzloc(i-1)))*(newtezhengt(i)-newtezhengt(i-1))/(tezhengt(i)-tezhengt(i-1))+newt(tzloc(i-1));
%     end 
% end

% f=f+0.1*rand(size(f));
CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=2;
Lambda = LSCurFit_V(P,k,N,newt,CList);

curve=LSMatrix_V(k,N,newt)*Lambda;
% err=vecnorm(P-curve,2,2);
% meanerr=mean(err);
% maxerr=max(err);
[errV,meanerrV,maxerrV,stderrV]=distanceerror(P,curve);




for i=2:length(tzloc)-1
    for j=0:1
        tempCList=CList;
        tempCList(ismember(CList(:,1),newtezhengt(i)),3)=j;
        tempLambda=LSCurFit_V(P,k,N,newt,tempCList);
        tempcurve=LSMatrix_V(k,N,newt)*tempLambda;
%         temperr=vecnorm(P-tempcurve,2,2);
%         tempmeanerr=mean(temperr);
        [temperrV,tempmeanerrV,tempmaxerrV,tempstderrV]=distanceerror(P,tempcurve);
        if  tempmeanerrV<meanerrV
            CList(ismember(CList(:,1),newtezhengt(i)),3)=j;
            Lambda=tempLambda;            
            meanerrV=tempmeanerrV;            
        end
    end
end
% LambdaS=Lambda;
% LambdaS(all(abs(Lambda)<=10^(-3),2),:)=0;
% LambdaS(all(LambdaS~=0,2),3)=find(all(abs(LambdaS)~=0,2));
% LambdaS(all(LambdaS==0,2),:)=[];
elapsedTimeV=toc;

% if plan<20
% load(['point',num2str(plan),'-200','.mat'])
% end
path='C:\CodeStore\matlab\vfitcurve\nosharpfit2D\';
Lambda(all(abs(Lambda)<=10^(-3),2),:)=0;
not0row=sum(any(Lambda,2));
%% 本方法最终重构结果
figure
curve=LSMatrix_V(k,N,newt)*Lambda;
% err=vecnorm(P-curve,2,2);
% maxerr=max(err);%最大误差
% meanerr=mean(err);%平均误差
% stderr=std(err);%标准差
% MSE=immse(P,curve);%均方误差
[errV,meanerrV,maxerrV,stderrV]=distanceerror(P,curve);
energyV=Lambda(:,1).^2+Lambda(:,2).^2;



plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10,'linewidth',1.5);hold on
% plot(curve(:,1),curve(:,2),'Color',[0 102 153]/255,'LineWidth',3)
VCompose(Lambda,k,N)
% plot(t,err,'linewidth',3)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% title('最终重构结果','fontsize', 15, 'fontname', '微软雅黑')
legend({'原始数据','重构曲线'},'fontsize', 15, 'fontname', '微软雅黑')
axis equal
axis off

% plot(P(tzloc,1),P(tzloc,2),'s','MarkerSize',10)
% plot(P(tzloc,1),P(tzloc,2),'s','MarkerSize',10)
duandian=LSMatrix_V(k,N,CList(:,1))*Lambda;
% plot(duandian(:,1),duandian(:,2),'s','MarkerSize',10)
saveas(gcf,[path,'plan',num2str(plan),'V','.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'V'])
print(gcf,'-dpng',[path,'plan',num2str(plan),'V'])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'V'])
%% 信号压缩结果
CList(end+1,:)=[0 1 0];
range2cut=1:32;

[~,slidx]=sort(abs(energyV),'descend');
% partlambda=Lambda(slidx(range2cut),:);
% partmat=LSMatrix_V(k,N,newt);
% partmat=partmat(:,slidx(range2cut));
% partcurve=partmat*partlambda;

partlambda=zeros(size(Lambda));
partlambda(slidx(range2cut),:)=Lambda(slidx(range2cut),:);
partlambda = LSCurFit_L2LC(partlambda,k,N,CList);
partcurve=LSMatrix_V(k,N,newt)*partlambda;
% errS=vecnorm(P-partcurve,2,2);
% maxerrS=max(errS);%最大误差
% meanerrS=mean(errS);%平均误差
% stderrS=std(errS);%标准差
[errS,meanerrS,maxerrS,stderrS]=distanceerror(P,partcurve);

figure
plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10,'linewidth',1.5);hold on
% plot(partcurve(:,1),partcurve(:,2),'Color',[0 102 153]/255,'LineWidth',3)
VCompose(partlambda,k,N)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend('原始数据','重构曲线','fontsize', 15, 'fontname', '微软雅黑')
axis equal
axis off
saveas(gcf,[path,'plan',num2str(plan),'VS',num2str(max(range2cut)),'.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'VS',num2str(max(range2cut))])
print(gcf,'-dpng',[path,'plan',num2str(plan),'VS',num2str(max(range2cut))])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'VS',num2str(max(range2cut))])


% %% 信号压缩结果2
% 
% frontN=4;
% front2cut=1:(k+1)*2^(frontN-1);
% partlambda2=Lambda(front2cut,:);
% fCList=zeros(2^(frontN-1)-1,3);
% fCList(:,1)=1/(2^(frontN-1)):1/(2^(frontN-1)):(2^(frontN-1)-1)/(2^(frontN-1));
% fCList(:,2)=fCList(:,1);
% fCList(:,3)=2;
% fCList(end,:)=[0 1 0];
% 
% partlambda2 = LSCurFit_L2LC(partlambda2,k,frontN,fCList);
% partmat2=LSMatrix_V(k,N,newt);
% partmat2=partmat2(:,front2cut);
% partcurve2=partmat2*partlambda2;
% % errS2=vecnorm(P-partcurve,2,2);
% % maxerrS=max(errS2);%最大误差
% % meanerrS=mean(errS2);%平均误差
% % stderrS=std(errS2);%标准差
% [errS2,meanerrS2,maxerrS2,stderrS2]=distanceerror(P,partcurve2);
% 
% figure
% plot(gpoint(:,1),gpoint(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
% plot(partcurve2(:,1),partcurve2(:,2),'Color',[0 102 153]/255,'LineWidth',3)
% %% PIA
% for i=1:runtime
% tic
% [ps,ix] = dpsimplify(P,fenbianlv);   %%0.0001
% Pp=[ps;gpoint(1,:)];
% k = 3;
% n = k + 1;
% pro = 1;
% NumItr =3;
% [PP_all,T] = PIA_CurApr(Pp(:,1:2),k,pro,NumItr);
% U = linspace(0,1,10*length(gpoint));                  % B样条曲线采样点对应参数值
% C1 = bspline_deboor(n,T,[real(PP_all(:,NumItr+1)),imag(PP_all(:,NumItr+1))],U);% 迭代NumItr步后曲线
% elapsedTimeP(i)=toc;
% end
% elapsedTimeP=mean(elapsedTimeP);
% figure
% plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10) 
% hold on
% plot(C1(:,1),C1(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% legend({'原始数据','重构曲线'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
% Chord = vecnorm(diff([gpoint;gpoint(1,:)],1,1),2,2);
% normt = [0;cumsum(Chord)/sum(Chord)];
% Cd = bspline_deboor(n,T,[real(PP_all(:,NumItr+1)),imag(PP_all(:,NumItr+1))],normt);
% [errP,meanerrP,maxerrP,stderrP]=distanceerror(P,Cd);
% axis equal
% axis off
% saveas(gcf,[path,'plan',num2str(plan),'PIA','.fig'])
% print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'PIA'])
% print(gcf,'-dpng',[path,'plan',num2str(plan),'PIA'])
% print(gcf,'-dsvg',[path,'plan',num2str(plan),'PIA'])
%% db4小波变换对比实验
figure
for i=1:runtime
tic
% nexttile
leveldb4=2;
cutdb4level=3;
[cxdb4,lxdb4] = wavedec(P(:,1),leveldb4,'db4');
rxdb4=cxdb4;
rxdb4(lxdb4(cutdb4level)+1:end) = 0;
iwxdb4 = waverec(rxdb4,lxdb4,'db4');

[cydb4,lydb4] = wavedec(P(:,2),leveldb4,'db4');
rydb4=cydb4;
rydb4(lxdb4(cutdb4level)+1:end) = 0;
iwydb4 = waverec(rydb4,lydb4,'db4');
elapsedTimedb4(i)=toc;
end
elapsedTimedb4=mean(elapsedTimedb4);
Pdb4=[iwxdb4 iwydb4];
Pdb4(end+1,:)=Pdb4(1,:);
energydb4=cxdb4.^2+cydb4.^2;
% errdb4=vecnorm(P-Pdb4);
% maxerrdb4=max(errdb4);%最大误差
% meanerrdb4=mean(errdb4);%平均误差
% stderrdb4=std(errdb4);%标准差
[errdb4,meanerrdb4,maxerrdb4,stderrdb4]=distanceerror(P,Pdb4);

plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10);hold on
plot(Pdb4(:,1),Pdb4(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)

set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend('原始数据','重构曲线','fontsize', 15, 'fontname', '微软雅黑')
axis equal
axis off
% title([num2str(ldb4(1)),'项db4小波变换'],'fontsize', 15, 'fontname', '微软雅黑')

saveas(gcf,[path,'plan',num2str(plan),'db4',num2str(sum(lxdb4(1:1+cutdb4level))),'.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'db4',num2str(sum(lxdb4(1:1+cutdb4level)))])
print(gcf,'-dpng',[path,'plan',num2str(plan),'db4',num2str(sum(lxdb4(1:1+cutdb4level)))])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'db4',num2str(sum(lxdb4(1:1+cutdb4level)))])
% %% db2小波变换对比实验
% figure
% leveldb2=3;
% cutdb2level=1;
% % nexttile
% [cxdb2,lxdb2] = wavedec(P(:,1),leveldb2,'db2');
% rxdb2=cxdb2;
% rxdb2(sum(lxdb2(1:1+cutdb2level))+1:end) = 0;
% iwxdb2 = waverec(rxdb2,lxdb2,'db2');
% 
% [cydb2,lydb2] = wavedec(P(:,2),leveldb2,'db2');
% rydb2=cydb2;
% rydb2(sum(lxdb2(1:1+cutdb2level))+1:end) = 0;
% iwydb2 = waverec(rydb2,lydb2,'db2');
% 
% Pdb2=[iwxdb2 iwydb2];
% 
% energydb2=iwxdb2.^2+iwydb2.^2;
% % errdb2=vecnorm(P-Pdb2);
% % maxerrdb2=max(errdb2);%最大误差
% % meanerrdb2=mean(errdb2);%平均误差
% % stderrdb2=std(errdb2);%标准差
% [errdb2,meanerrdb2,maxerrdb2,stderrdb2]=distanceerror(P,Pdb2);
% 
% plot(gpoint(:,1),gpoint(:,2),'.','Color',[255 102 102]/255,'MarkerSize',10);hold on
% plot(Pdb2(:,1),Pdb2(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)
% 
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑','looseInset',get(gca,'TightInset'))
% legend('原始数据点','特征点','fontsize', 15, 'fontname', '微软雅黑')
% axis equal
% axis off
% % title([num2str(ldb2(1)),'项db2小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
% 
% saveas(gcf,[path,'plan',num2str(plan),'db2',num2str(sum(lxdb2(1:1+cutdb2level))),'.fig'])
% print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'db2',num2str(sum(lxdb2(1:1+cutdb2level)))])
% print(gcf,'-dpng',[path,'plan',num2str(plan),'db2',num2str(sum(lxdb2(1:1+cutdb2level)))])
% print(gcf,'-dsvg',[path,'plan',num2str(plan),'db2',num2str(sum(lxdb2(1:1+cutdb2level)))])


%% coif

figure
for i=1:runtime
tic
% nexttile
levelcoif3=2;
cutcoif3level=3;
[cxcoif3,lxcoif3] = wavedec(P(:,1),levelcoif3,'coif3');
rxcoif3=cxcoif3;
rxcoif3(lxcoif3(cutcoif3level)+1:end) = 0;
iwxcoif3 = waverec(rxcoif3,lxcoif3,'coif3');

[cycoif3,lycoif3] = wavedec(P(:,2),levelcoif3,'coif3');
rycoif3=cycoif3;
rycoif3(lxcoif3(cutcoif3level)+1:end) = 0;
iwycoif3 = waverec(rycoif3,lycoif3,'coif3');
sum(lxcoif3(1:1+cutcoif3level))
Pcoif3=[iwxcoif3 iwycoif3];
elapsedTimecoif3(i)=toc;
end
elapsedTimecoif3=mean(elapsedTimecoif3);
energycoif3=cxcoif3.^2+cycoif3.^2;
% errcoif3=vecnorm(P-Pcoif3);
% maxerrcoif3=max(errcoif3);%最大误差
% meanerrcoif3=mean(errcoif3);%平均误差
% stderrcoif3=std(errcoif3);%标准差
[errcoif3,meanerrcoif3,maxerrcoif3,stderrcoif3]=distanceerror(P,Pcoif3);

Pcoif3(end+1,:)=Pcoif3(1,:);
plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10);hold on
plot(Pcoif3(:,1),Pcoif3(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)

set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend('原始数据','重构曲线','location','northwest','fontsize', 15, 'fontname', '微软雅黑')
axis equal
axis off
% title([num2str(lcoif3(1)),'项coif3小波变换'],'fontsize', 15, 'fontname', '微软雅黑')

saveas(gcf,[path,'plan',num2str(plan),'coif3',num2str(sum(lxcoif3(1:1+cutcoif3level))),'.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'coif3',num2str(sum(lxcoif3(1:1+cutcoif3level)))])
print(gcf,'-dpng',[path,'plan',num2str(plan),'coif3',num2str(sum(lxcoif3(1:1+cutcoif3level)))])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'coif3',num2str(sum(lxcoif3(1:1+cutcoif3level)))])
%% sym5

figure
for i=1:runtime
tic
% nexttile
levelsym5=2;
cutsym5level=3;
[cxsym5,lxsym5] = wavedec(P(:,1),levelsym5,'sym5');
rxsym5=cxsym5;
rxsym5(lxsym5(cutsym5level)+1:end) = 0;
iwxsym5 = waverec(rxsym5,lxsym5,'sym5');

[cysym5,lysym5] = wavedec(P(:,2),levelsym5,'sym5');
rysym5=cysym5;
rysym5(lxsym5(cutsym5level)+1:end) = 0;
iwysym5 = waverec(rysym5,lysym5,'sym5');
sum(lxsym5(1:1+cutsym5level))
Psym5=[iwxsym5 iwysym5];
elapsedTimesym5(i)=toc;
end
elapsedTimesym5=mean(elapsedTimesym5);
energysym5=cxsym5.^2+cysym5.^2;
% errsym5=vecnorm(P-Psym5);
% maxerrsym5=max(errsym5);%最大误差
% meanerrsym5=mean(errsym5);%平均误差
% stderrsym5=std(errsym5);%标准差
[errsym5,meanerrsym5,maxerrsym5,stderrsym5]=distanceerror(P,Psym5);

Psym5(end+1,:)=Psym5(1,:);
plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10);hold on
plot(Psym5(:,1),Psym5(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)

set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend('原始数据','重构曲线','fontsize', 15, 'fontname', '微软雅黑')
axis equal
axis off
% title([num2str(lsym5(1)),'项sym5小波变换'],'fontsize', 15, 'fontname', '微软雅黑')

saveas(gcf,[path,'plan',num2str(plan),'sym5',num2str(sum(lxsym5(1:1+cutsym5level))),'.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'sym5',num2str(sum(lxsym5(1:1+cutsym5level)))])
print(gcf,'-dpng',[path,'plan',num2str(plan),'sym5',num2str(sum(lxsym5(1:1+cutsym5level)))])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'sym5',num2str(sum(lxsym5(1:1+cutsym5level)))])

%% Fourier变换对比试验
figure
for i=1:runtime
tic
% nexttile
nd=189;
z = frdescp(P);
P_FFT = ifrdescp(z,nd);
elapsedTimef(i)=toc;
end
elapsedTimef=mean(elapsedTimef);
energyf=real(z).^2+imag(z).^2;
% errf=abs(P-P_FFT);
% maxerrf=max(errf);%最大误差
% meanerrf=mean(errf);%平均误差
% stderrf=std(errf);%标准差
[errf,meanerrf,maxerrf,stderrf]=distanceerror(P,P_FFT);

plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10);hold on
plot(P_FFT(:,1),P_FFT(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)


set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend('原始数据','重构曲线','fontsize', 15, 'fontname', '微软雅黑')
axis equal
axis off

% title([num2str(nd),'项Fourier变换'],'fontsize', 15, 'fontname', '微软雅黑')

saveas(gcf,[path,'plan',num2str(plan),'f',num2str(nd),'.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'f',num2str(nd)])
print(gcf,'-dpng',[path,'plan',num2str(plan),'f',num2str(nd)])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'f',num2str(nd)])


%% 能量集中度
range=9:113;
ceV=cumsum(sort(energyV(range),'descend'));
normveV=ceV/ceV(end);
cedb4=cumsum(sort(energydb4(range),'descend'));
normvedb4=cedb4/cedb4(end);
% cedb2=cumsum(sort(energydb2(5:max(range2cut)),'descend'));
% normvedb2=cedb2/cedb2(end);
ceF=cumsum(sort(energyf(range),'descend'));
normveF=ceF/ceF(end);
cecoif3=cumsum(sort(energycoif3(range),'descend'));
normvecoif3=cecoif3/cecoif3(end);
cesym5=cumsum(sort(energysym5(range),'descend'));
normvesym5=cesym5/cesym5(end);
figure
hold on
plot(normveV,'LineWidth',3)
plot(normveF,'LineWidth',3)
plot(normvedb4,'LineWidth',3)
plot(normvecoif3,'LineWidth',3)
plot(normvesym5,'LineWidth',3)
% plot(normvedb2,'LineWidth',3)
% plot(normenergystl(1:256),'LineWidth',1.5)
legend({'本文方法','Fourier','Db4','Coif3','Sym5'},'location','best','fontsize', 15, 'fontname', '微软雅黑')

ylim([0 1.2])
xlim([0 length(range)])
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')

saveas(gcf,[path,'plan',num2str(plan),'NL','.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'NL'])
print(gcf,'-dpng',[path,'plan',num2str(plan),'NL'])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'NL'])


%% 3D全局误差图
figure
% load(['errB',num2str(plan),'.mat'])

x=linspace(1,length(errV),length(errV));
y=ones(size(errV));
hold on
box on
plot3(x,1*y,errV,'LineWidth',1.1);
plot3(x,2*y,errf,'LineWidth',1.1);
plot3(x,3*y,errdb4,'LineWidth',1.1);
plot3(x,4*y,errcoif3,'LineWidth',1.1);
plot3(x,5*y,errsym5,'LineWidth',1.1);
% plot3(x,2*y,errP,'LineWidth',1.1);
% plot3(x,3*y,errf,'LineWidth',1.1);
% plot3(x,4*y,errdb4,'LineWidth',1.1);
% plot3(x,5*y,errcoif3,'LineWidth',1.1);
% plot3(x,6*y,errsym5,'LineWidth',1.1);
% plot3(x,7*y,wuchaB,'LineWidth',1.1);
for i=1:length(tzloc)
    plot3(x(tzloc(i))*y,linspace(0,7,length(errV)),0*y,'k','linewidth',1.5);
end

view(11,17)
xlim([1 length(errV)]);
ylim([0.8 5.2]);
% zlabel('global error')
zlabel('全局误差')
% set(gca,'ytick',[1: 4],'yticklabel',{'V-system','Genetic algorithm','Fourier','Db4 wavelet'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% set(gca,'ytick',[1: 3],'yticklabel',{'V-system','Fourier','Db4 wavelet'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
set(gca,'ytick',[1:5],'yticklabel',{'本文方法','Fourier','Db4','Coif3','Sym5'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')

% set(gca,'ytick',[1:7],'yticklabel',{'本章方法','PIA','Fourier','Db4','Coif3','Sym5','GA'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
saveas(gcf,[path,'plan',num2str(plan),'ER3D','.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'ER3D'])
print(gcf,'-dpng',[path,'plan',num2str(plan),'ER3D'])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'ER3D'])


function [phi,K]=qulv(gpoint,i,num)
r=i+1;
l=i-1;
if i==1
    l=num-1;
elseif i==num
    r=2;
end
v1=gpoint(i,:)-gpoint(l,:);
v2=gpoint(r,:)-gpoint(i,:);
v3=gpoint(r,:)-gpoint(l,:);
phi=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
K=2*sin(phi)/norm(v3);
end

function y=xiaoshu2dec(x,N)
y=zeros(1,N+1);
tempx=x;
for i=1:N+1
    
    y(i)=floor(tempx*2);
    tempx=tempx*2-y(i);
end
if y(end)==1
    y(N)=y(N)+1;
    y(end)=[];
end
for i=fliplr(2:N)
    if y(i)==2
        y(i)=0;
        y(i-1)=y(i-1)+1;
    end
end
% if y(1)==2
%     for i=1:N
%         y(i)=floor(x*2);
%         x=x*2-y(i);
%     end
% end
end

function y=dec2xiaoshu(x)
for i=1:length(x)
    temp(i)=x(i)*2^(-i);
end
y=sum(temp);
end
function [kappa,norm_k] = PJcurvature(gpoint,i)
    x = gpoint(i-1:i+1,1);
    y = gpoint(i-1:i+1,2);
    t_a = norm([x(2)-x(1),y(2)-y(1)]);
    t_b = norm([x(3)-x(2),y(3)-y(2)]);
    
    M =[[1, -t_a, t_a^2];
        [1, 0,    0    ];
        [1,  t_b, t_b^2]];

    a = M\x;
    b = M\y;

    kappa  = 2.*abs(a(3)*b(2)-b(3)*a(2)) / (a(2)^2.+b(2)^2.)^(1.5);
    norm_k =  [b(2),-a(2)]/sqrt(a(2)^2.+b(2)^2.);
end
