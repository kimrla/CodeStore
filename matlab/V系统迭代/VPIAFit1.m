close all %不局部叠加，用局部基函数重新求解方程
clear
example=5;
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
        R = loadsvg('clef_lspia.svg', 0.01,0 );
        C = unique(R{1,1},'rows','stable');
        
        P=curvspace(C,200);
    case 14
        % 实验1
        x=0:0.005:1-0.005;
        f=90./(1+exp(-100*(x-0.5)));
        %         M=length(x);%采样点个数
        %         N=32;%基函数个数
        P=[x*100;f]';
    case 15
        % 实验2
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-5))+(10*x-5).^5/500;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        P=[x*100;f]';
    case 12
        % 实验2 变体3
        x=0:0.02:1-0.02;
        f=100./exp(abs(10*x-5))+(10*x-5).^5/30;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        P=[x*100;f]';
    case 11
        % 实验2 变体2
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-5));
        P=[x*100;f]';
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
    case 13
        % 实验2 变体2
        x=0:0.0002:1-0.0002;
        f=100./exp(abs((10*x-2.5).*(10*x-7.5)));
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        P=[x*100;f]';
    case 10
        % 实验2 变体1
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-6))+(10*x-6).^5/500;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        P=[x*100;f]';
    case 3
        %
        x=0:0.005:1;
        f=1./(0.01+(x-0.3).^2).*(x<0.5)+1./(0.015+(x-0.65).^2).*(x>=0.5);
        %         M=length(x);
        %         N=32;%基函数个数
        P=[x*100;f]';
    case 4
        load fire500.mat
        P=gpoint;
        P=[P;P(1,:)];
    case 5
        load yezi600.mat
        P=gpoint;
        P=[P;P(1,:)];
    case 6
        load G-200.mat
        P=gpoint;
        P=[P;P(1,:)];
    case 7
        load hudie2fjy420.mat
        P=gpoint;
        P=[P;P(1,:)];
        
end

pointnum=length(P);
scatter(P(:,1),P(:,2));
axis equal
t=canshuhua(P);
% %%
% 
% for i=1:pointnum
%     [phi(i),K(i)]=qulv(P,i,pointnum);
% end
% q=zeros(pointnum,1);
% for i=1:pointnum
%     r=i+1;
%     l=i-1;
%     if i==1
%         l=pointnum-1;
%     elseif i==pointnum
%         r=2;
%     end
%     if K(i)>K(r) && K(i)>K(l) && phi(i)>pi/6
%         %         && phi(i)>pi/6        &&K(i)>3*mean(K)
%         q(i)=1;
%     end
% end
% tezhengt=t(q==1);
% Nleast=ceil(log2(length(tezhengt)-1)+1);
% for i=1:length(nprate)
%     newnprate(i)=dec2xiaoshu(xiaoshu2dec(nprate(i),N));
% end
% while length(newnprate)-length(unique(newnprate))
%     N=N+1;
%     for i=1:length(nprate)
%         newnprate(i)=dec2xiaoshu(xiaoshu2dec(nprate(i),N));
%     end
%     
% end
       

%%
% load fire500newt.mat
% load fire500newtzt.mat
% t=newt;
% %%
% k=3;
% N=5;
% V=LSMatrix_V(k,N,t);
% Lambda0=V\P;
% scatter(P(:,1),P(:,2));hold on
% VCompose(Lambda0,k,N)
% figure
% V(t>=0.5 & t<=1,9:64)=0;
% Lambda1=V\P;
% scatter(P(:,1),P(:,2));hold on
% VCompose(Lambda1,k,N)





%%
k=3;
N=1;%初始V组数


tol1 = 1; % 允许的最大误差
tol2 = 1e-6; % 要求的最小误差变化
err1 = tol1 + 1;%误差初始值
err2 = tol2 + 1;%误差变化初始值
gen=0;%当前迭代数
Gmax=200;%最大迭代数
maxout=tol1+1;%局部片段误差初始值
tpiece=[0,1];%初始参数区间划分
tidx=[1,pointnum];%初始区间划分索引

% tt=linspace(0,1,100*pointnum)';
% Pp=interp1(t,P,tt,'linear');
Pp=P;
P0=Pp;
Lambda1=[];
Nmax=floor(log2(pointnum/(k+1)))+1;
while(err1 > tol1 && err2 > tol2 &&err2>0 && gen<Gmax&&N<Nmax)
    
    fidx=find(maxout<=tol1);%找出不需要加层的区间索引
    
    N=N+1;
    gen=gen+1;
    vmat=LSMatrix_V(k,N,t);    
    
    A = LSMatrix_V(k,N,t);
    clear CList
    CList=zeros(2^(N-1)-1,3);
    CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
    CList(:,2)=CList(:,1);
    CList(:,3)=2;
%     CList(ismember(CList(:,1),newtezhengt),3)=0;
    
    if N==2
%         Lambda0=zeros((k+1)*2^(N-1),2);
        temA=A;
        Lambda1=pinv(temA)*Pp;
    else
        
        temB=A(:,(k+1)*2^(N-2)+1:end);%取最高层基函数
        
        for i=1:length(fidx)    
            temB(t>=tpiece(fidx(i))&t<=tpiece(fidx(i)+1),:)=0;
            
        end
            temA=[temA,temB];
            Lambda1=pinv(temA)*Pp;
        Lambda1 = LSCurFit_L2LC(Lambda1,k,N,CList);
    end
        
    
        
        Lambda=pinv(A)*P0;%标准V系统拟合解
        delta0=vecnorm(P-vmat*Lambda,2,2);%标准解的各点误差
        maxerr0=max(delta0);
        meanerr0=mean(delta0);
    
    deltaP=P- vmat*Lambda1;%拟合曲线和原始点的差
    delta = vecnorm(deltaP,2,2);%每个点的误差
    
%     Pp=Pp-A*Lambda1;
    tmp=max(delta);
    meanerr=mean(delta);%平均误差
    err2=abs(err1-tmp);
    err1=tmp;
    tpiece=linspace(0,1,2^(N-1)+1);%参数分割区间
    tidx=[];
    for i=1:2^(N-1)+1
        tidx(end+1)=find(t>=tpiece(i),1);
    end
    
    m=max(diff(tidx));
    clear out
    out(1:m,1:length(tidx)-1)=nan;
    
    for i=1:length(tidx)-2
        out(1:tidx(i+1)-tidx(i),i)=delta(tidx(i):tidx(i+1)-1);
    end
        i=i+1;
        out(1:tidx(i+1)-tidx(i)+1,i)=delta(tidx(i):tidx(i+1));
    [maxout,maxoutidx]=max(out,[],'omitnan');
    
% %     
    figure,
    plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
    VCompose(Lambda1,k,N)
%     figure
%     plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
%     VCompose(Lambda,k,N)
%     VCompose(Lambdatrapz,k,N)
end
% Lambda1(all(abs(Lambda1)<=10^(-3),2),:)=0;
% Lambda1(all(Lambda1~=0,2),3)=find(all(Lambda1~=0,2));
% Lambda1(all(Lambda1==0,2),:)=[];
% delta=0.0001;
% CList=zeros(2^(N-1)-1,3);
% CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
% CList(:,2)=CList(:,1);
% CList(:,3)=2;
Lambda = LSCurFit_V(P,k,N,t,CList);
% Lambda1 = LSCurFit_L2LC(Lambda1,k,N,CList);
fdelta=vecnorm(P-vmat*Lambda1,2,2);
fmaxerr=max(fdelta);
fmeanerr=mean(fdelta);
fdelta0=vecnorm(P-vmat*Lambda,2,2);
fmaxerr0=max(fdelta0);
fmeanerr0=mean(fdelta0);
% Lambdatrapz=LSCurFit_Vtrapz(P,k,N,t,CList);
% Lambda=LSMatrix_V(k,N,t)\P;
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
figure,
plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
VCompose(Lambda1,k,N)
figure,
plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
VCompose(Lambda,k,N)
% VCompose(Lambdatrapz,k,N)
function [CList,Lambdaupdate]=VPIA(j,k,N,t,P,CList)

lt=(j-1)/2^(N-2);
rt=j/2^(N-2);

CList(end+1,:)=[(lt+rt)/2,(lt+rt)/2,2];
utidx=find((lt<=t)&(t<=rt));
tupdate=t(utidx);
Pupdate=P(utidx,:);
Lambdaupdate=LSCurFit_VPIA(Pupdate,k,N,tupdate,CList,j);
%     Lambdaupdate=LSMatrix_V(k,N,tupdate)\Pupdate;
CList=CList(end,:);
end

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
