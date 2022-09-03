clear %局部加层解方程时扩展区间
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
        R = loadsvg('clef_lspia.svg', 0.01,0 );
        P = unique(R{1,1},'rows','stable');
        pointnum=length(P);
end
scatter(P(:,1),P(:,2));
axis equal
t=canshuhua(P);

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
    
    fidx=find(maxout>tol1);%找出需要加层的区间索引
    
    N=N+1;
    gen=gen+1;
    vmat=LSMatrix_V(k,N,t);    
    
    A = LSMatrix_V(k,N,t);
    if N==2
%         Lambda0=zeros((k+1)*2^(N-1),2);
        temA=A;
        Lambda0=pinv(temA)*Pp;
    else
        Lambda0=zeros((k+1)*2^(N-2),2);
%         temA=A(:,(k+1)*2^(N-2)+1:end);
    
        for i=1:length(fidx)  
            
            temt=2^(N-2)*(t(t>=tpiece(fidx(i))&t<=tpiece(fidx(i)+1))-tpiece(fidx(i)));%把加层区间的参数映射到[0,1]
            temA=sqrt(2^(N-2))*GenerateValue(k,temt);%本层基函数
            
            Lambda0(i:2^(N-2):i+k*2^(N-2),:)=pinv(temA)*Pp(t>=tpiece(fidx(i))&t<=tpiece(fidx(i)+1),:);
        end
    end
        Lambda1=[Lambda1;Lambda0];
    
        
        Lambda=pinv(A)*P0;
    
    deltaP=P- vmat*Lambda1;%拟合曲线和原始点的差
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
    figure
    plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
    VCompose(Lambda,k,N)
%     VCompose(Lambdatrapz,k,N)
end
Lambda1(all(Lambda1~=0,2),3)=find(all(Lambda1~=0,2));
Lambda1(all(Lambda1==0,2),:)=[];
% delta=0.0001;
CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=2;
Lambda = LSCurFit_V(P,k,N,t,CList);
Lambdatrapz=LSCurFit_Vtrapz(P,k,N,t,CList);
% Lambda=LSMatrix_V(k,N,t)\P;
Pp=LSMatrix_V(k,N,t)*Lambda;
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
VCompose(Lambda,k,N)
VCompose(Lambdatrapz,k,N)
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

    function G = GenerateValue(k,x)
    % 生成元函数在x上的取值
    % 注意：当x=0.5时的取值（左右极限之算术平均值）
        thed = 10^(-10);
        x1_idx = find(x>=0 & (0.5-x)>thed);
        x2_idx = find((x-0.5)>thed & x<=1);
        x12_idx = [x1_idx;x2_idx];
        x1 = x(x1_idx);
        x2 = x(x2_idx);
        NumX = length(x);
        G = zeros(NumX,k+1);
        GValue_mid = [1,k+1];  % 生成元在x=0.5处取值
        switch k
            case 1
                G(x12_idx,:) = [sqrt(3)*[-4 * x1 + 1;  4 * x2 - 3],...
                             [-6 * x1 + 1; -6 * x2 + 5]];
                GValue_mid = [-sqrt(3),0];         
            case 2
                G(x12_idx,:) = [sqrt(5)*[16 * x1.^2 - 10 * x1 + 1; -16 * x2.^2 + 22 * x2 -  7],...
                     sqrt(3)*[30 * x1.^2 - 14 * x1 + 1;  30 * x2.^2 - 46 * x2 + 17],...
                             [40 * x1.^2 - 16 * x1 + 1; -40 * x2.^2 + 64 * x2 - 25]];
                GValue_mid = [0,sqrt(27)/2,0];
            case 3
                G(x12_idx,:) = [sqrt(7)*[ -64 * x1.^3 +  66 * x1.^2 - 18 * x1 + 1;   64 * x2.^3 - 126 * x2.^2 +  78 * x2 - 15],...
                     sqrt(5)*[-140 * x1.^3 + 114 * x1.^2 - 24 * x1 + 1; -140 * x2.^3 + 306 * x2.^2 - 216 * x2 + 49],...
                     sqrt(3)*[-224 * x1.^3 + 156 * x1.^2 - 28 * x1 + 1;  224 * x2.^3 - 516 * x2.^2 + 388 * x2 - 95],...
                             [-280 * x1.^3 + 180 * x1.^2 - 30 * x1 + 1; -280 * x2.^3 + 660 * x2.^2 - 510 * x2 + 129]];
                 GValue_mid = [sqrt(7)/2,0,-2*sqrt(3),0];
        end
        
        x05_idx = find(abs(x-0.5)<=thed);
        if ~isempty(x05_idx)
            G(x05_idx,:) =  repmat(GValue_mid,[length(x05_idx),1]);
        end
            
%         if NumX ~= length(x12_idx)
%             if ~isempty(find(x==0))
%                G(find(x==0),:) =  GValue_bp(1,:);
%             end
%             if ~isempty(find(x==0.5))
%                G(find(x==0.5),:) =  GValue_bp;
%             end
%             if ~isempty(find(x==1))
%                G(find(x==1),:) =  GValue_bp(3,:);
%             end
%         end
        
    end