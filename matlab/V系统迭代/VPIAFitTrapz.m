close %局部加层用积分
clear
example=2;
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
        pointnum=length(P);
end

scatter(P(:,1),P(:,2));hold on
axis equal


t=canshuhua(P);

k=3;
N=1;%初始V组数


tol1 = 4; % 允许的最大误差
tol2 = 1e-6; % 要求的最小误差变化
err1 = tol1 + 1;%误差初始值
err2 = tol2 + 1;%误差变化初始值
gen=0;%当前迭代数
Gmax=200;%最大迭代数
maxout=tol1+1;%局部片段误差初始值
tpiece=[0,1];%初始参数区间划分
tidx=[1,pointnum];%初始区间划分索引

Lambda1=[];
Nmax=floor(log2(pointnum/(k+1)));
% tt=linspace(0,1,10*pointnum)';
tt=linspace(0,1,2^(Nmax+8)+1)';
tt=unique([t;tt]);
Pp=interp1(t,P,tt,'linear');
plot(Pp(:,1),Pp(:,2),'.')
P0=Pp;

while(err1 > tol1 && err2 > tol2 && gen<Gmax&&N<Nmax)
    
    fidx=find(maxout>tol1);%找出需要加层的区间索引
    
    N=N+1;
    gen=gen+1;
    vmat=LSMatrix_V(k,N,t);    
%     clear CList
%     CList=zeros(2^(N-1)-1,3);
%     CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
%     CList(:,2)=CList(:,1);
%     CList(:,3)=2;
    A = LSMatrix_V(k,N,tt);
    if N==2
        Lambda0=zeros((k+1)*2^(N-1),2);
        temA=A;
    else
        Lambda0=zeros((k+1)*2^(N-2),2);
        temA=A(:,(k+1)*2^(N-2)+1:end);
    end
    for i=1:length(fidx)
        temP=zeros(length(Pp),2);
        if ismember(fidx(i)+1,fidx)
            temP(tt>=tpiece(fidx(i))&tt<tpiece(fidx(i)+1),:)=Pp(tt>=tpiece(fidx(i))&tt<tpiece(fidx(i)+1),:);
        else
            temP(tt>=tpiece(fidx(i))&tt<=tpiece(fidx(i)+1),:)=Pp(tt>=tpiece(fidx(i))&tt<=tpiece(fidx(i)+1),:);
        end
        for j=1:size(P,2)
            l=1:size(temA,2);
            Lambda0(l,j)=trapz(tt,temP(:,j).*temA(:,l))'+Lambda0(l,j);%用数值积分形式分层求lambda1
        end
        
    end
        Lambda1=[Lambda1;Lambda0];
%         figure,
%         plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
%         VCompose(Lambda1,k,N)
%         Lambda1 = LSCurFit_L2LC(Lambda1,k,N,CList);
        
    
        for j=1:size(P,2)
            l=1:size(A,2);
            Lambda(l,j)=trapz(tt,P0(:,j).*A(:,l));%用数值积分形式求lambda
        end
    
        Lambda2=vmat\P;
    
    deltaP=P- vmat*Lambda1;%拟合曲线和原始数据点的差向量
    delta = vecnorm(deltaP,2,2);%每个点的误差
    
    Pp=Pp-A*Lambda1;
    tmp=max(delta);
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
    maxout=max(out,[],'omitnan');
% %     
    figure,
    plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
    VCompose(Lambda1,k,N)
    figure,
    plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
    VCompose(Lambda,k,N)
    figure,
    plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
    VCompose(Lambda2,k,N)
%     VCompose(Lambdatrapz,k,N)
end

CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=2;
Lambda1 = LSCurFit_L2LC(Lambda1,k,N,CList);
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

figure,
plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
VCompose(Lambda1,k,N)
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
