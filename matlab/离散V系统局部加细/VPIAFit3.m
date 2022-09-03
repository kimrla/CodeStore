close all%局部加层用正交V矩阵求解
clear all
plan=11;
switch plan
    case 1
        fid=fopen('C:\Users\J\桌面\coord_seligFmt\coord_seligFmt\ag35.dat','r');
        dataCell = textscan(fid,'%f %f','headerlines',1);
        
        P = cell2mat(dataCell);
        % if P(end,:)~=P(1,:)
        %     P(end+1,:)=P(1,:);
        % end
        % P=P+1;
        pointnum=250;
        csl = [0;cumsum(vecnorm(diff(P),2,2))]; %累加长度
        L=csl(end);
        
        t=csl/L;
        tt=linspace(0,1,pointnum)';
        
        Pp=interp1(t,P,tt,'linear');
        
        
        clear P
        P=Pp;
        clear Pp
    case 2
        R = loadsvg('clef_lspia.svg', 0.1,0 );
        P = unique(R{1,1},'rows','stable');
        %         pointnum=length(P);
    case 3
        %
        x=0:0.005:1;
        f=1./(0.01+(x-0.3).^2).*(x<0.5)+1./(0.015+(x-0.65).^2).*(x>=0.5);
        %         M=length(x);
        %         N=32;%基函数个数
        P=[x*100;f]';
    case 4
        load star3fjy360.mat
    case 5
        load hainan1987.mat
    case 6
        load ecg101_2_512.mat
    case 7
        load hudie2fjy420.mat
    case 8
        load G-200.mat
        P=gpoint;
        P(end+1,:)=P(1,:);
    case 9
        load niao21000.mat
        P=gpoint;
        P(end+1,:)=P(1,:);
    case 10
        load 苏州市-1024.mat
        P=P{1,1};
        P=[P(58:end,:);P(1:57,:)];
%         P=[P(580:end,:);P(1:57,:)];
%         [ps,ix] = dpsimplify(P,10);
%         P=P(1/2*length(P)+5:end,:);
        fenbianlv=0.05;
    case 11
        load 巴彦淖尔市-1024.mat
        P=P{1,1};
        P=[P(903:end,:);P(1:902,:)];
        fenbianlv=0.1;
    case 12
        load dgx2_256.mat
    case 13
        load jianlong-1024.mat
        P=P{1,1};
        fenbianlv=0.02;
    case 14
        load zhuazi-1024.mat
        P=P{1,1};
        fenbianlv=0.02;
    case 15 
        load cat1-1024.mat
        P=P{1,1};
        fenbianlv=0.08;
    case 16
        load baolong-1024.mat
        P=P{1,1};
        fenbianlv=0.2;
    case 17
        load sanjiaolong-1024.mat
        P=P{1,1};
        fenbianlv=0.02;
    case 18
        load bat-1024.mat
        P=P{1,1};
        fenbianlv=0.1;
    case 19
        load yezishu-1024.mat
        P=P{1,1};
        fenbianlv=0.06;
    case 20
        load 白城市-1024.mat
        P=P{1,1};
        P=[P(99:end,:);P(1:98,:)];
        fenbianlv=0.08;
end
runtime=1;
scatter(P(:,1),P(:,2));hold on
axis equal

tic
t=canshuhua(P);
k=3;
Nmax=9;

num=(k+1)*2^(Nmax-1);
tt=linspace(1,2*num-1,num)'/(2*num);%对原始数据插值取点需要的均匀参数
PP=interp1(t,P,tt,'linear');%对原始数据插值取点
% PP=P;
A=LSMatrix_V(k,Nmax,tt);
% O=Schmidt_orthogonalization(A);
matrixname=['V',num2str(k),'_',num2str(num),'.mat'];
load(matrixname)
O(abs(O)<1e-15)=0;


% save(matrixname,'O')
% % A=orth(A);
% E=O*O';

N=1;%初始V组数

tol1 = 1e-1; % 允许的最大误差
tol2 = 1e-6; % 要求的最小误差变化
tol3=1e-3;%允许的最大平均误差

err1 = tol1 + 1;%误差初始值
err2 = tol2 + 1;%误差变化初始值
meanerr=tol3+1;%平均误差初始值
meanout=meanerr;
gen=0;%当前迭代数

maxout=tol1+1;%局部片段误差初始值
tpiece=[0,1];%初始参数区间划分
tidx=[1,num];%初始区间划分索引

Lambda1=[];
figure
plot(PP(:,1),PP(:,2),'.','Color','k','MarkerSize',10)
P0=PP;
axis equal
axis off
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')

while(err1 > tol1  &&N<Nmax)
    fidx=find(maxout>tol1);%找出需要加层的区间索引
    % while(meanerr > tol3 && gen<Gmax&&N<Nmax)
    %     fidx=find(meanout>tol3);%找出需要加层的区间索引
    N=N+1;
    gen=gen+1;
    vmat=LSMatrix_V(k,N,t);
    %     clear CList
    %     CList=zeros(2^(N-1)-1,3);
    %     CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
    %     CList(:,2)=CList(:,1);
    %     CList(:,3)=2;
    %     A = LSMatrix_V(k,N,tt);
    if N==2
        
        temA=O(:,1:(k+1)*2);
        Lambda0=temA'*PP;
    else
        Lambda0=zeros((k+1)*2^(N-2),2);
        temA=O(:,(k+1)*2^(N-2)+1:(k+1)*2^(N-1));
        
        for i=1:length(fidx)
            temP=PP(tt>tpiece(fidx(i))&tt<tpiece(fidx(i)+1),:);
            Lambda0(fidx(i):2^(N-2):fidx(i)+k*2^(N-2),:)=temA(tt>tpiece(fidx(i))&tt<tpiece(fidx(i)+1),fidx(i):2^(N-2):fidx(i)+k*2^(N-2))'*temP;
        end
    end
    Lambda1=[Lambda1;Lambda0];
    %         figure,
    %         plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
    %         VCompose(Lambda1,k,N)
    %         Lambda1 = LSCurFit_L2LC(Lambda1,k,N,CList);
    
    
    
            Lambda3=O(:,1:(k+1)*2^(N-1))'*P0-Lambda1;
            
    %         Lambda2=vmat\P;
    PC=O(:,1:(k+1)*2^(N-1))*Lambda1;
    deltaP=P0- PC;%拟合曲线和原始数据点的差向量
    delta = vecnorm(deltaP,2,2);%每个点的误差
    meanerr=mean(delta);
    PP=deltaP;
    tmp=max(delta);
    err2=abs(err1-tmp);
    err1=tmp;
    tpiece=linspace(0,1,2^(N-1)+1);%参数分割区间
    tidx=[];
    for i=1:2^(N-1)
        tidx(end+1)=find(tt>=tpiece(i),1);
    end
    i=i+1;
    tidx(end+1)=num;
    
    m=max(diff(tidx));
    clear out
    out(1:m,1:length(tidx)-1)=nan;
    
    for i=1:length(tidx)-2
        out(1:tidx(i+1)-tidx(i),i)=delta(tidx(i):tidx(i+1)-1);
    end
    i=i+1;
    out(1:tidx(i+1)-tidx(i)+1,i)=delta(tidx(i):tidx(i+1));
    maxout=max(out,[],'omitnan');
    meanout=mean(out,'omitnan');
    % %
%     figure,
%     plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
%     %     VCompose(Lambda1,k,N)
%     
%     plot(PC(:,1),PC(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
%     
%     %     LSVCompose(Lambda1,k,Nmax,O)
%     
%     %     figure,
%     %     plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
%     %     VCompose(Lambda2,k,N)
%     %     VCompose(Lambdatrapz,k,N)
%     %     fig1name=['OVPIA',num2str(N),num2str(1/tol1)];
%     %     savefig(fig1name)
%     %     saveas(gcf,fig1name,'png')
%     %     saveas(gcf,fig1name,'emf')
%     set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
%     title(['第',num2str(N),'代变换结果'],'fontsize', 15, 'fontname', '微软雅黑')
%     axis equal
%     axis off
end
elapsedTimeV=toc;
% CList=zeros(2^(N-1)-1,3);
% CList(:,1)=(1:(2^(N-1)-1))/(2^(N-1));
% CList(:,2)=CList(:,1);
% CList(:,3)=0;
% Lambda1 = LSCurFit_L2LC(Lambda1,k,N,CList);
% Lambda1(all(Lambda1~=0,2),3)=find(all(Lambda1~=0,2));
% Lambda1(all(Lambda1==0,2),:)=[];

% delta=0.0001;
% CList=zeros(2^(N-1)-1,3);
% CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
% CList(:,2)=CList(:,1);
% CList(:,3)=2;
% Lambda1=LSCurFit_L2LC(Lambda1,k,N,CList);
% figure,
% plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
% VCompose(Lambda1,k,N)
%
% Lambda = LSCurFit_V(P,k,N,t,CList);
% Lambdatrapz=LSCurFit_Vtrapz(P,k,N,t,CList);
% % Lambda=LSMatrix_V(k,N,t)\P;
% Pp=LSMatrix_V(k,N,t)*Lambda;
% wucha=vecnorm((P-Pp),2,2);
% % Lambdaupdate=zeros((k+1)*2^(N-2),1);
% N=N+1;
% NumSeg = 2^(N-2);
% for j=1:NumSeg
%     [CList1,Lambdaupdate(j:NumSeg:size(Lambda),:)]=VPIA(j,k,N,t,P-Pp,CList);
% end
% Lambda=[Lambda;Lambdaupdate];
% Aupdate=LSMatrix_V(k,N,t);
% Aupdate=Aupdate(:,size(Aupdate,2)/2+1:end);

% PpP=Aupdate*Lambdaupdate;
path='C:\CodeStore\matlab\V系统迭代\pic\';
Lambda1(all(abs(Lambda1)<=20*10^(-2),2),:)=0;
% Lambda(all(Lambda~=0,2),3)=find(all(Lambda~=0,2));
% Lambda(all(Lambda==0,2),:)=[];
%% 本方法重构结果
figure,
plot(P(:,1),P(:,2),'.','Color','r','MarkerSize',10);hold on
curve=O*Lambda1;

plot(curve(:,1),curve(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)
% notzeronum=length(find(all(Lambda1~=0,2)));
[errV,meanerrV,maxerrV,stderrV] = distanceerror(P0,curve);
curve(end+1,:)=P(1,:);
energyV=Lambda1(:,1).^2+Lambda1(:,2).^2;
MSE=mean(errV.^2);
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend({'原始数据','拟合曲线'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
% title('最终重构结果','fontsize', 15, 'fontname', '微软雅黑')
not0row=sum(any(Lambda1,2));
    axis equal
    axis off
saveas(gcf,[path,'plan',num2str(plan),'V','.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'V'])
print(gcf,'-dpng',[path,'plan',num2str(plan),'V'])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'V'])
% [error,ave,maxda,scha] = distanceerror(P,cPP);
%     fig1name=['OVPIA',num2str(N),num2str(1/tol1)];
%     savefig(fig1name)
%     saveas(gcf,fig1name,'png')
%     saveas(gcf,fig1name,'emf')



% %% %和V系统做对比
% lamdaV=LSCurFit_V(P0,k,N,tt,CList);
% figure
% plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
% % VCompose(lamdaV,k,N)
% cPV=LSMatrix_V(k,N,tt)*lamdaV;
% % cPV=LSMatrix_V(k,N,t)*lamdaV;
% [errorV,aveV,maxdaV,schaV] = distanceerror(P0,cPV);
% % [errorV,aveV,maxdaV,schaV] = distanceerror(P,cPV);
% plot(cPV(:,1),cPV(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
% %     fig1name=['V',num2str(N),num2str(1/tol1)];
% %     savefig(fig1name)
% %     saveas(gcf,fig1name,'png')
% %     saveas(gcf,fig1name,'emf')
% %% 和正交V做对比
% figure
% LAMDAOV=O(:,1:(k+1)*2^(N-1))'*P0;
% LAMDAOV=LSCurFit_L2LC(LAMDAOV,k,N,CList);
% CURVEOV=O(:,1:(k+1)*2^(N-1))*LAMDAOV;
% CPOV=interp1(tt,CURVEOV,t,'makima');
% [errorOV,aveOV,maxdaOV,schaOV] = distanceerror(P0,CURVEOV);
% % [errorOV,aveOV,maxdaOV,schaOV] = distanceerror(P,CPOV);
% plot(P0(:,1),P0(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
% plot(CPOV(:,1),CPOV(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
% %     fig1name=['OV',num2str(Nmax),num2str(1/tol1)];
% %     savefig(fig1name)
% %     saveas(gcf,fig1name,'png')
% %     saveas(gcf,fig1name,'emf')
%% 信号压缩结果

range2cut=1:32;
[~,slidx]=sort(abs(energyV),'descend');
partlambda=zeros(size(Lambda1));
partlambda(slidx(range2cut),:)=Lambda1(slidx(range2cut),:);
partcurve=O*partlambda;

[errS,meanerrS,maxerrS,stderrS]=distanceerror(P,partcurve);
partcurve(end+1,:)=P0(1,:);
figure
plot(P0(:,1),P0(:,2),'.','Color','r','MarkerSize',10,'linewidth',1.5);hold on
plot(partcurve(:,1),partcurve(:,2),'Color',[0 102 153]/255,'LineWidth',2)

set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend('原始数据','拟合曲线','location','best','fontsize', 15, 'fontname', '微软雅黑')
axis equal
axis off
saveas(gcf,[path,'plan',num2str(plan),'VS',num2str(max(range2cut)),'.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'VS',num2str(max(range2cut))])
print(gcf,'-dpng',[path,'plan',num2str(plan),'VS',num2str(max(range2cut))])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'VS',num2str(max(range2cut))])
%% 省略高频后的结果
% figure
% lowhz=1/4*length(Lambda1);
% curveJ=O(:,1:lowhz)*Lambda1(1:lowhz,:);
% plot(P0(:,1),P0(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
% plot(curveJ(:,1),curveJ(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
% legend({'原始数据','重构曲线'},'fontsize', 15, 'fontname', '微软雅黑')
% CPJ=interp1(tt,curveJ,t,'makima');hold on
% [errorJ,aveJ,maxdaJ,schaJ] = distanceerror(P0,curveJ);
% MSEJ=mean(errorJ.^2);
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% title(['前',num2str(lowhz),'项V基函数重构结果'],'fontsize', 15, 'fontname', '微软雅黑')
%     axis equal
%     axis off

% [errorJ,aveJ,maxdaJ,schaJ] = distanceerror(P,CPJ);
%     fig1name=['OVPIAJ',num2str(N),num2str(1/tol1)];
%     savefig(fig1name)
%     saveas(gcf,fig1name,'png')
%     saveas(gcf,fig1name,'emf')
% 
% figure
% MatV=LSMatrix_V(k,N,tt);
% % MatV=LSMatrix_V(k,N,t);
% cPVJ=MatV(:,1:lowhz)*lamdaV(1:lowhz,:);
% [errorVJ,aveVJ,maxdaVJ,schaVJ] = distanceerror(P,cPVJ);
% plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
% plot(cPVJ(:,1),cPVJ(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
% 
% %     fig1name=['VJ',num2str(N),num2str(1/tol1)];
% %     savefig(fig1name)
% %     saveas(gcf,fig1name,'png')
% %     saveas(gcf,fig1name,'emf')
% figure
% curveOVJ=O(:,1:size(LAMDAOV(1:lowhz,:),1))*LAMDAOV(1:lowhz,:);
% plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
% plot(curveOVJ(:,1),curveOVJ(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
% CPOVJ=interp1(tt,curveOVJ,t,'makima');hold on
% [errorOVJ,aveOVJ,maxdaOVJ,schaOVJ] = distanceerror(P0,curveOVJ);
% [errorOVJ,aveOVJ,maxdaOVJ,schaOVJ] = distanceerror(P,CPOVJ);
%     fig1name=['OVJ',num2str(N),num2str(1/tol1)];
%     savefig(fig1name)
%     saveas(gcf,fig1name,'png')
%     saveas(gcf,fig1name,'emf')

%% PIA
% fenbianlv=2.4;
for i=1:runtime
tic
[ps,ix] = dpsimplify(P0,fenbianlv);   %%0.0001
Pp=[ps;P0(1,:)];
k = 3;
n = k + 1;
pro = 1;
NumItr =5;
[PP_all,T] = PIA_CurApr(Pp(:,1:2),k,pro,NumItr);
U = linspace(0,1,10*length(P0));                  % B样条曲线采样点对应参数值
C1 = bspline_deboor(n,T,[real(PP_all(:,NumItr+1)),imag(PP_all(:,NumItr+1))],U);% 迭代NumItr步后曲线
elapsedTimeP(i)=toc;
end
elapsedTimeP=mean(elapsedTimeP);
figure
plot(P0(:,1),P0(:,2),'.','Color','r','MarkerSize',10) 
hold on
plot(C1(:,1),C1(:,2),'Color',[0 102 153]/255,'LineWidth',2)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend({'原始数据','拟合曲线'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
Chord = vecnorm(diff([P0;P0(1,:)],1,1),2,2);
normt = [0;cumsum(Chord)/sum(Chord)];
Cd = bspline_deboor(n,T,[real(PP_all(:,NumItr+1)),imag(PP_all(:,NumItr+1))],normt);
[errP,meanerrP,maxerrP,stderrP]=distanceerror(P0,Cd);
axis equal
axis off
saveas(gcf,[path,'plan',num2str(plan),'PIA',num2str(length(ps)),'.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'PIA',num2str(length(ps))])
print(gcf,'-dpng',[path,'plan',num2str(plan),'PIA',num2str(length(ps))])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'PIA',num2str(length(ps))])
%% db4小波变换对比实验
figure
for i=1:runtime
tic
leveldb4=2;
cutdb4level=3;
[cxdb4,lxdb4] = wavedec(P0(:,1),leveldb4,'db4');
cxdb4(lxdb4(cutdb4level)+1:end) = 0;
iwxdb4 = waverec(cxdb4,lxdb4,'db4');
[cydb4,lydb4] = wavedec(P0(:,2),leveldb4,'db4');
cydb4(lydb4(cutdb4level)+1:end) = 0;
iwydb4 = waverec(cydb4,lydb4,'db4');
P_DWTdb4 = [iwxdb4,iwydb4];
elapsedTimedb4(i)=toc;
end
elapsedTimedb4=mean(elapsedTimedb4);
energydb4=cxdb4.^2+cydb4.^2;
[errdb4,meanerrdb4,maxerrdb4,stderrdb4] = distanceerror(P0,P_DWTdb4);
P_DWTdb4(end+1,:)=P0(1,:);
MSEdb4=mean(errdb4.^2);
plot(P0(:,1),P0(:,2),'.','Color','r','MarkerSize',10);hold on
plot(P_DWTdb4(:,1),P_DWTdb4(:,2),'Color',[0 102 153]/255,'LineWidth',2)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend({'原始数据','拟合曲线'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
% title([num2str(lxdb4(3)),'项db4小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
    axis equal
    axis off
saveas(gcf,[path,'plan',num2str(plan),'DB4',num2str(lxdb4(cutdb4level)),'.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'DB4',num2str(lxdb4(cutdb4level))])
print(gcf,'-dpng',[path,'plan',num2str(plan),'DB4',num2str(lxdb4(cutdb4level))])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'DB4',num2str(lxdb4(cutdb4level))])
%% db2小波变换对比实验
% figure
% [cxdb2,lxdb2] = wavedec(P0(:,1),2,'db2');
% cxdb2(lxdb2(1)+1:end) = 0;
% iwxdb2 = waverec(cxdb2,lxdb2,'db2');
% [cydb2,lxdb2] = wavedec(P0(:,2),2,'db2');
% cydb2(lxdb2(1)+1:end) = 0;
% iwydb2 = waverec(cydb2,lxdb2,'db2');
% P_DWTdb2 = [iwxdb2,iwydb2];
% [errorwdb2,avewdb2,maxdawdb2,schawdb2] = distanceerror(P,P_DWTdb2);
% MSEdb2=mean(errorwdb2.^2);
% plot(P0(:,1),P0(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
% plot(P_DWTdb2(:,1),P_DWTdb2(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% title([num2str(lxdb2(1)),'项db4小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
%     axis equal
%     axis off
%% Fourier变换对比试验
figure
for i=1:runtime
clear z
tic
% nexttile
nd=500;
z = frdescp(P0);
P_FFT = ifrdescp(z,nd);
elapsedTimef(i)=toc;
end
elapsedTimef=mean(elapsedTimef);
energyf=real(z).^2+imag(z).^2;
% errf=abs(P-P_FFT);
% maxerrf=max(errf);%最大误差
% meanerrf=mean(errf);%平均误差
% stderrf=std(errf);%标准差
[errf,meanerrf,maxerrf,stderrf]=distanceerror(P0,P_FFT);

plot(P0(:,1),P0(:,2),'.','Color','r','MarkerSize',10);hold on
plot(P_FFT(:,1),P_FFT(:,2),'Color',[0 102 153]/255,'LineWidth',2)


set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
legend('原始数据','拟合曲线','location','best','fontsize', 15, 'fontname', '微软雅黑')
axis equal
axis off

% title([num2str(nd),'项Fourier变换'],'fontsize', 15, 'fontname', '微软雅黑')

saveas(gcf,[path,'plan',num2str(plan),'f',num2str(nd),'.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'f',num2str(nd)])
print(gcf,'-dpng',[path,'plan',num2str(plan),'f',num2str(nd)])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'f',num2str(nd)])

% % coif
% 
% figure
% for i=1:runtime
% tic
% nexttile
% levelcoif3=2;
% cutcoif3level=3;
% [cxcoif3,lxcoif3] = wavedec(P0(:,1),levelcoif3,'coif3');
% rxcoif3=cxcoif3;
% rxcoif3(lxcoif3(cutcoif3level)+1:end) = 0;
% iwxcoif3 = waverec(rxcoif3,lxcoif3,'coif3');
% 
% [cycoif3,lycoif3] = wavedec(P0(:,2),levelcoif3,'coif3');
% rycoif3=cycoif3;
% rycoif3(lxcoif3(cutcoif3level)+1:end) = 0;
% iwycoif3 = waverec(rycoif3,lycoif3,'coif3');
% 
% Pcoif3=[iwxcoif3 iwycoif3];
% elapsedTimecoif3(i)=toc;
% end
% elapsedTimecoif3=mean(elapsedTimecoif3);
% energycoif3=cxcoif3.^2+cycoif3.^2;
% errcoif3=vecnorm(P-Pcoif3);
% maxerrcoif3=max(errcoif3);%最大误差
% meanerrcoif3=mean(errcoif3);%平均误差
% stderrcoif3=std(errcoif3);%标准差
% [errcoif3,meanerrcoif3,maxerrcoif3,stderrcoif3]=distanceerror(P0,Pcoif3);
% 
% Pcoif3(end+1,:)=Pcoif3(1,:);
% plot(P0(:,1),P0(:,2),'.','Color','r','MarkerSize',10);hold on
% plot(Pcoif3(:,1),Pcoif3(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)
% 
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% legend('原始数据','拟合曲线','location','best','fontsize', 15, 'fontname', '微软雅黑')
% axis equal
% axis off
% title([num2str(lcoif3(1)),'项coif3小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
% 
% saveas(gcf,[path,'plan',num2str(plan),'coif3',num2str(sum(lxcoif3(1:1+cutcoif3level))),'.fig'])
% print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'coif3',num2str(lxcoif3(cutcoif3level))])
% print(gcf,'-dpng',[path,'plan',num2str(plan),'coif3',num2str(lxcoif3(cutcoif3level))])
% print(gcf,'-dsvg',[path,'plan',num2str(plan),'coif3',num2str(lxcoif3(cutcoif3level))])
% % sym5
% figure
% for i=1:runtime
% tic
% 
% nexttile
% levelsym5=2;
% cutsym5level=3;
% [cxsym5,lxsym5] = wavedec(P0(:,1),levelsym5,'sym5');
% rxsym5=cxsym5;
% rxsym5(lxsym5(cutsym5level)+1:end) = 0;
% iwxsym5 = waverec(rxsym5,lxsym5,'sym5');
% 
% [cysym5,lysym5] = wavedec(P0(:,2),levelsym5,'sym5');
% rysym5=cysym5;
% rysym5(lxsym5(cutsym5level)+1:end) = 0;
% iwysym5 = waverec(rysym5,lysym5,'sym5');
% 
% Psym5=[iwxsym5 iwysym5];
% elapsedTimesym5(i)=toc;
% end
% elapsedTimesym5=mean(elapsedTimesym5);
% energysym5=cxsym5.^2+cysym5.^2;
% errsym5=vecnorm(P-Psym5);
% maxerrsym5=max(errsym5);%最大误差
% meanerrsym5=mean(errsym5);%平均误差
% stderrsym5=std(errsym5);%标准差
% [errsym5,meanerrsym5,maxerrsym5,stderrsym5]=distanceerror(P0,Psym5);
% 
% Psym5(end+1,:)=Psym5(1,:);
% plot(P0(:,1),P0(:,2),'.','Color','r','MarkerSize',10);hold on
% plot(Psym5(:,1),Psym5(:,2),'Color',[0 102 153]/255,'LineWidth',1.1)
% 
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% legend('原始数据','拟合曲线','location','best','fontsize', 15, 'fontname', '微软雅黑')
% axis equal
% axis off
% title([num2str(lsym5(1)),'项sym5小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
% 
% saveas(gcf,[path,'plan',num2str(plan),'sym5',num2str(lxsym5(cutsym5level)),'.fig'])
% print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'sym5',num2str(lxsym5(cutsym5level))])
% print(gcf,'-dpng',[path,'plan',num2str(plan),'sym5',num2str(lxsym5(cutsym5level))])
% print(gcf,'-dsvg',[path,'plan',num2str(plan),'sym5',num2str(lxsym5(cutsym5level))])


%% 能量集中度
range=9:300;
ceV=cumsum(sort(energyV(range),'descend'));
normveV=ceV/ceV(end);
cedb4=cumsum(sort(energydb4(range),'descend'));
normvedb4=cedb4/cedb4(end);
% cedb2=cumsum(sort(energydb2(5:max(range2cut)),'descend'));
% normvedb2=cedb2/cedb2(end);
ceF=cumsum(sort(energyf(range),'descend'));
normveF=ceF/ceF(end);
% cecoif3=cumsum(sort(energycoif3(range),'descend'));
% normvecoif3=cecoif3/cecoif3(end);
% cesym5=cumsum(sort(energysym5(range),'descend'));
% normvesym5=cesym5/cesym5(end);
figure
hold on
plot(normveV,'LineWidth',3)
plot(normveF,'LineWidth',3)
plot(normvedb4,'LineWidth',3)
% plot(normvecoif3,'LineWidth',3)
% plot(normvesym5,'LineWidth',3)
% plot(normvedb2,'LineWidth',3)
% plot(normenergystl(1:256),'LineWidth',1.5)
% legend({'本章方法','Fourier','Db4','Coif3','Sym5'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
legend({'本章方法','Fourier','Db4'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
xlim([0 max(range)])
ylim([0 1.2])
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')

saveas(gcf,[path,'plan',num2str(plan),'NL','.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'NL'])
print(gcf,'-dpng',[path,'plan',num2str(plan),'NL'])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'NL'])

%% 3D全局误差图
figure
% load(['errB',num2str(plan),'.mat'])
% wuchaB(end)=[];
x=linspace(1,length(errV),length(errV));
y=ones(size(errV));
hold on
box on
plot3(x,1*y,errV,'LineWidth',1.1);
plot3(x,2*y,errf,'LineWidth',1.1);
plot3(x,3*y,errdb4,'LineWidth',1.1);
plot3(x,4*y,errP,'LineWidth',1.1);
% plot3(x,4*y,errcoif3,'LineWidth',1.1);
% plot3(x,5*y,errsym5,'LineWidth',1.1);
% plot3(x,2*y,errP,'LineWidth',1.1);
% plot3(x,3*y,errf,'LineWidth',1.1);
% plot3(x,4*y,errdb4,'LineWidth',1.1);
% plot3(x,5*y,errcoif3,'LineWidth',1.1);
% plot3(x,6*y,errsym5,'LineWidth',1.1);
% plot3(x,7*y,wuchaB,'LineWidth',1.1);
% for i=1:length(tzloc)
%     plot3(x(tzloc(i))*y,linspace(0,5,length(errV)),0*y,'k','linewidth',1.5);
% end

view(11,17)
xlim([1 length(errV)]);
ylim([0.8 4.2]);
% zlabel('global error')
zlabel('全局误差')
% set(gca,'ytick',[1: 4],'yticklabel',{'V-system','Genetic algorithm','Fourier','Db4 wavelet'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
set(gca,'ytick',[1: 5],'yticklabel',{'本章方法','Fourier','Db4','PIA'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% set(gca,'ytick',[1:7],'yticklabel',{'本章方法','PIA','Fourier','Db4','Coif3','Sym5','GA'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
saveas(gcf,[path,'plan',num2str(plan),'ER3D','.fig'])
print(gcf,'-depsc','-tiff',[path,'plan',num2str(plan),'ER3D'])
print(gcf,'-dpng',[path,'plan',num2str(plan),'ER3D'])
print(gcf,'-dsvg',[path,'plan',num2str(plan),'ER3D'])
