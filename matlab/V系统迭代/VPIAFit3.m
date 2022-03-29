close all%局部加层用正交V矩阵求解
clear
example=10;
switch example
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
        [ps,ix] = dpsimplify(P,10);
    case 11
        load 巴彦淖尔市-1024.mat
        P=P{1,1};
    case 12
        load dgx2_256.mat
        
end

scatter(P(:,1),P(:,2));hold on
axis equal


t=canshuhua(P);
k=3;
Nmax=9;

num=(k+1)*2^(Nmax-1);
tt=linspace(1,2*num-1,num)'/(2*num);%对原始数据插值取点需要的均匀参数
% PP=interp1(t,P,tt,'linear');%对原始数据插值取点
% PP=P;
A=LSMatrix_V(k,Nmax,tt);
O=Schmidt_orthogonalization(A);
O(abs(O)<1e-15)=0;
matrixname=['V',num2str(k),'_',num2str(num),'.mat'];
% load(matrixname)
% save(matrixname,'O')
% % A=orth(A);
% E=O*O';

N=1;%初始V组数

tol1 = 1e-3; % 允许的最大误差
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
plot(PP(:,1),PP(:,2),'.')
P0=PP;

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
            Lambda0(i:2^(N-2):i+k*2^(N-2),:)=temA(tt>tpiece(fidx(i))&tt<tpiece(fidx(i)+1),i:2^(N-2):i+k*2^(N-2))'*temP;
        end
    end
    Lambda1=[Lambda1;Lambda0];
    %         figure,
    %         plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
    %         VCompose(Lambda1,k,N)
    %         Lambda1 = LSCurFit_L2LC(Lambda1,k,N,CList);
    
    
    
    %         Lambda3=O(:,1:(k+1)*2^(N-1))'*P0;
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
    figure,
    plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
    %     VCompose(Lambda1,k,N)
    
    plot(PC(:,1),PC(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
    
    %     LSVCompose(Lambda1,k,Nmax,O)
    
    %     figure,
    %     plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
    %     VCompose(Lambda2,k,N)
    %     VCompose(Lambdatrapz,k,N)
    %     fig1name=['OVPIA',num2str(N),num2str(1/tol1)];
    %     savefig(fig1name)
    %     saveas(gcf,fig1name,'png')
    %     saveas(gcf,fig1name,'emf')
    set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
    title(['第',num2str(N),'代变换结果'],'fontsize', 15, 'fontname', '微软雅黑')
end

CList=zeros(2^(N-1)-1,3);
CList(:,1)=(1:(2^(N-1)-1))/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=0;
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
%% 本方法重构结果
figure,
plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
curve=O*Lambda1;
plot(curve(:,1),curve(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
notzeronum=length(find(all(Lambda1~=0,2)));
[error,ave,maxda,scha] = distanceerror(P0,curve);
MSE=mean(error.^2);
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('最终重构结果','fontsize', 15, 'fontname', '微软雅黑')
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
%% 本方法将基函数系数排列取有限个最大值重构结果
figure
[~,slidx]=sort(sum(Lambda1.^2,2),'descend');
firstnum=256;
sortedlam=zeros(length(Lambda1),2);
sortedlam(slidx(1:firstnum),:)=Lambda1(slidx(1:firstnum),:);
curvestl=O*sortedlam;
[errorstl,avestl,maxdastl,schastl] = distanceerror(P0,curvestl);
MSEstl=mean(errorstl.^2);
plot(P0(:,1),P0(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(curvestl(:,1),curvestl(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title(['系数最大',num2str(firstnum),'项V基函数重构结果'],'fontsize', 15, 'fontname', '微软雅黑')

%% 省略高频后的结果
figure
lowhz=1/4*length(Lambda1);
curveJ=O(:,1:lowhz)*Lambda1(1:lowhz,:);
plot(P0(:,1),P0(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(curveJ(:,1),curveJ(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
CPJ=interp1(tt,curveJ,t,'makima');hold on
[errorJ,aveJ,maxdaJ,schaJ] = distanceerror(P0,curveJ);
MSEJ=mean(errorJ.^2);
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title(['前',num2str(lowhz),'项V基函数重构结果'],'fontsize', 15, 'fontname', '微软雅黑')

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
%% db4小波变换对比实验
figure
[cxdb4,lxdb4] = wavedec(P0(:,1),2,'db4');
cxdb4(lxdb4(1)+1:end) = 0;
iwxdb4 = waverec(cxdb4,lxdb4,'db4');
[cydb4,lydb4] = wavedec(P0(:,2),2,'db4');
cydb4(lydb4(1)+1:end) = 0;
iwydb4 = waverec(cydb4,lydb4,'db4');
P_DWTdb4 = [iwxdb4,iwydb4];
[errorwdb4,avewdb4,maxdawdb4,schawdb4] = distanceerror(P,P_DWTdb4);
MSEdb4=mean(errorwdb4.^2);
plot(P0(:,1),P0(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(P_DWTdb4(:,1),P_DWTdb4(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(lxdb4(1)),'项db4小波变换'],'fontsize', 15, 'fontname', '微软雅黑')
%% db2小波变换对比实验
figure
[cxdb2,lxdb2] = wavedec(P0(:,1),2,'db2');
cxdb2(lxdb2(1)+1:end) = 0;
iwxdb2 = waverec(cxdb2,lxdb2,'db2');
[cydb2,lxdb2] = wavedec(P0(:,2),2,'db2');
cydb2(lxdb2(1)+1:end) = 0;
iwydb2 = waverec(cydb2,lxdb2,'db2');
P_DWTdb2 = [iwxdb2,iwydb2];
[errorwdb2,avewdb2,maxdawdb2,schawdb2] = distanceerror(P,P_DWTdb2);
MSEdb2=mean(errorwdb2.^2);
plot(P0(:,1),P0(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(P_DWTdb2(:,1),P_DWTdb2(:,2),'Color',[0 102 153]/255,'LineWidth',1.5)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(lxdb2(1)),'项db4小波变换'],'fontsize', 15, 'fontname', '微软雅黑')

%% Fourier变换对比试验
figure
nd=256;
zx = frdescp1d(P0(:,1));
P_FFTx = ifrdescp1d(zx,nd);
zy = frdescp1d(P0(:,2));
P_FFTy = ifrdescp1d(zy,nd);
P_FFT=[P_FFTx,P_FFTy];
[errorf,avef,maxdaf,schaf] = distanceerror(P,P_FFT);
MSEf=mean(errorf.^2);
plot(P0(:,1),P0(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
plot(P_FFTx,P_FFTy,'Color',[0 102 153]/255,'LineWidth',1.5)
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title([num2str(nd),'项Fourier变换'],'fontsize', 15, 'fontname', '微软雅黑')

%% 累积归一化能量图
figure
energyv=cumsum(sort(vecnorm(Lambda1(4:260,:),2,2),'descend'));
normenergyv=energyv/energyv(end);
% energyv=cumsum(sort(abs(Lambda1(1:256,1))+abs(Lambda1(1:256,2)),'descend'));
% normenergyv=energyv/energyv(end);
energyf=cumsum(sort(vecnorm([abs(zx(1:256)),abs(zy(1:256))],2,2),'descend'));
normenergyf=energyf/energyf(end);
energydb4=cumsum(sort(vecnorm([cxdb4(1:256),cydb4(1:256)],2,2),'descend'));
nomenergydb4=energydb4/energydb4(end);
energydb2=cumsum(sort(vecnorm([cxdb2(1:256),cydb2(1:256)],2,2),'descend'));
nomenergydb2=energydb2/energydb2(end);
% energystl=[0;cumsum(sort(sortedlam.^2,'descend'))];
% normenergystl=energystl/energystl(end);

hold on
plot(normenergyv,'LineWidth',1.5)
plot(normenergyf,'LineWidth',1.5)
plot(nomenergydb4,'LineWidth',1.5)
plot(nomenergydb2,'LineWidth',1.5)
% plot(normenergystl(1:256),'LineWidth',1.5)
legend({'V系统正交变换','Fourier变换','db4小波变换','db2小波变换'},'location','best','fontsize', 15, 'fontname', '微软雅黑')
xlim([0 256])
ylim([0 1.2])
set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
title('累积归一化频谱能量图','fontsize', 15, 'fontname', '微软雅黑')
