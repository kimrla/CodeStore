clear
example=13;
switch example
    case 1
        % 实验1
        x=0:0.005:1-0.005;
        f=90./(1+exp(-100*(x-0.5)));
        %         M=length(x);%采样点个数
        %         N=32;%基函数个数
        xmax=1;
    case 2
        % 实验2
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-5))+(10*x-5).^5/500;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        xmax=1;
    case 12
        % 实验2 变体3
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-5))+(10*x-5).^5/30;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        xmax=1;
    case 11
        % 实验2 变体2
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-5));
        xmax=1;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
    case 13
        % 实验2 变体2
        x=0:0.0002:1-0.0002;
        f=100./exp(abs((10*x-2.5).*(10*x-7.5)));
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        xmax=1;
    case 10
        % 实验2 变体1
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-6))+(10*x-6).^5/500;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        xmax=1;
    case 3
        %
        x=0:0.005:1-0.005;
        f=1./(0.01+(x-0.3).^2).*(x<0.5)+1./(0.015+(x-0.65).^2).*(x>=0.5);
        %         M=length(x);
        %         N=32;%基函数个数
        xmax=1;
    case 4
        x=0:0.02:10-0.02;
        y=@(x) 100./exp(abs(x-5))+(x-5).^5/500;
        wide=1;
        for i=1:length(x)
            f(i)=1/wide*integral(y,x(i)-wide/2,x(i)+wide/2);
        end
        xmax=10;
    case 5
        x=-2:0.04:2-0.04;
        f=sin(x)+2*exp(-30*x.^2);
        x=x+2;
        xmax=4;
    case 6
        x=-2:0.04:2-0.04;
        f=sin(2*x)+2*exp(-16*x.^2)+2;
        x=x+2;
        xmax=4;
    case 7
        x=0:0.005:1-0.005;
        f=4*x.^2.*(3-4*x).*(0<=x&x<0.5)+(4/3*x.*(4*x.^2-10*x+7)-3/2).*(0.5<=x&x<0.75)+(16/3*x.*(x-1).^2).*(0.75<=x&x<=1);
        xmax=1;
    case 8
        x=0:0.005:1-0.005;
        f=2*sin(4*pi*x)-6*(abs(x-0.4)).^0.3-0.5*sign(0.7-x);
        xmax=1;
    case 9
        x=0:0.005:1-0.005;
        f=sin(4*x-2)+2*exp(-30*(4*x-2).^2);
        xmax=1;
end

% f_=f+normrnd(0,1,1,M);

% k=1;%次数
% % t=canshuhua(M,d);
t=x/xmax;
t=t';
f=f';
% d=[t;f_]';
% [V,lamda]=vxishu(k,N,M,t,d);
% [V1,lamda1]=vdxishu(k,N,M,t,d);
% plot(t,f,'*')
% hold on
% drawV(N, k, lamda);
% % figure
% % for i=1:N
% %     subplot(N,1,i),plot(x,V(:,i))
P=[x',f];
% % end
k = 3;
N = 2;
% A = LSMatrix_V(k,N,t);

plan=1;


switch plan
    
    
    case 1
        % % Show
        % InfoV1Bas = BaseGene_V1(N);   % 1次V系统基函数信息(非离散采样)
        % [VRInfo,NumSeg] = VReconstruction_Polyline([X Y],InfoV1Bas);
        % figure,
        % plot(t,f,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
        % for j = 1 : NumSeg
        %     plot(VRInfo(j,5:6),VRInfo(j,7:8),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
        % end
        % % axis equal
        % % axis on
        num=length(P);
        for i=1:num
            [phiP(i),KP(i)]=qulv(P,i,num);
        end
        [qljdz,qljdd]=findpeaks(KP);
        qljdt=t(qljdd);
        
        for n=1:(k+1)*2^(N-1)
            A = LSMatrix_V(k,N,t);
            %         n=1;
            A=A(:,1:n);
            Lambda=A\P;
            X=Lambda(:,1);
            Y=Lambda(:,2);
            V = LSMatrix_V(k,N,linspace(0,1,1000)');
            V=V(:,1:n);
            C1=V*Lambda;
            Pp=A*Lambda;
            if n==1
                C1(:,1)=linspace(0,1,1000)';
                Pp(:,1)=x';
            end
            
            wucha=vecnorm((Pp-P),2,2);
            
            tz=0.5;
            for i=1:num
                
                [phiPp(i),KPp(i)]=qulv(Pp,i,num);
            end
            phiwc=phiPp./phiP;
            Kwc=KPp./KP;
            tzPhiwc(n)=phiwc(x==tz);
            tzKwc(n)=Kwc(x==tz);
            pj(n)=mean(wucha);
            tzwucha(n)=wucha(x==tz);
            % V = LSMatrix_V(k,N,linspace(0,1,501)');
            % F = V * X;
            % figure,
            % plot(t,f,'.r');hold on
            % plot(linspace(0,1,501)',F,'b');
            %         figure,
            subplot(2^(N-1),k+1,n),
            plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',10);hold on
            plot(C1(:,1),C1(:,2),'Color',[0 102 153]/255,'LineWidth',2)
            title(['添加前' num2str(n) '个基函数拟合效果'])
            
            %         axis equal
            %         axis off
            xlim([0 1])
            xlabel('x')
            ylabel('f(x)')
            box off
        end
        
        [rgx,rgy]=jianruijida(P,k,t);
        %         legend('原始数据','拟合曲线')
        for i=1:size(qljdt,1)
            for h=1:size(rgx,1)
                if qljdt(i,1)>=rgx(h,1) && qljdt(i,1)<=rgx(h,2)
                    qljdt(i,2)=h;                    
                end
            end
            for j=1:size(rgy,1)
                if qljdt(i,1)>=rgy(j,1) && qljdt(i,1)<=rgy(j,2)
                    qljdt(i,3)=j;
                end
            end
        end
        
        pjbhl=[0,-diff(pj)];
        figure
        plot(pj,'Color',[0 102 153]/255,'LineWidth',2)
        xlabel('V-系统基函数序号')
        ylabel('平均误差变化')
        figure
        plot(tzPhiwc,'Color',[0 102 153]/255,'LineWidth',2)
        xlabel('V-系统基函数序号')
        ylabel('特征点处拟合曲线与原始数据夹角比值')
        box off
    case 2
        % 约束方程
        %         [DR,DL] = VContinuityInfo1(N);
        
        % CList = [1/8,1/8,0;
        %             0.25,0.25,0;
        %             3/8,3/8,0;
        %             0.5,0.5,0;
        %             5/8,5/8,0;
        %             0.75,0.75,0;
        %             7/8,7/8,0;];
        CList=zeros(2^(N-1)-1,3);
        CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
        CList(:,2)=CList(:,1);
        CList(:,3)=0;
        %         tz=0.5;
        %         CList(ismember(CList(:,1),tz),:)=[];
        %         %       CList(16,:)=[];
        %         p = sum(CList(:,3)+1);  % 约束方程数目
        %         SegNum = 2^(N-1); % 分段数
        %         VNum = 2^N;       % 基函数数目
        %         knots = linspace(0,1,SegNum+1);% 节点向量
        %         C = zeros(p,VNum);    % 约束矩阵
        %         csidx = 0;
        %         for h = 1 : length(CList)           % 给每个约束建立相应的方程组
        %             x1 = find(knots==CList(h,1));   % 左节点
        %             x2 = find(knots==CList(h,2));   % 右节点
        %             for c = 0 : CList(h,3)          % 从C^0约束至C^CList(h,3)约束，每个约束建立一个方程
        %                 csidx = csidx + 1;
        %                 C(csidx,:) = [DR(:,x1,c+1) - DL(:,x2,c+1)]';
        %             end
        %         end
        %
        %         M = [2*A'*A,C';
        %             C,zeros(p)];
        %         d = zeros(p,1);
        %         bx = [2*A'*t;d];
        %         by = [2*A'*f;d];
        %
        %         X = M\bx;
        %         Y = M\by;
        %
        %         X = X(1:VNum);
        %         Y = Y(1:VNum);
        
        Lambda = LSCurFit_V(P,k,N,t,CList);
        X=Lambda(:,1);
        Y=Lambda(:,2);
        
        %% Show2
        figure,
        plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
        VCompose(Lambda,k,N)
        %         axis equal
        %         axis off
        xlim([0 1])
        xlabel('x')
        ylabel('f(x)')
        box off
    case 3
        X = A\t;
        Y=A\f;
        
        %         lastN=abs(Y((k+1)*2^(N-2)+1:(k+1)*2^(N-1)));
        for i=1:k+1
            lastN(:,i)=abs(Y(((k+1)*2^(N-2)+1+2^(N-2)*(i-1):(k+1)*2^(N-2)+2^(N-2)*i)));
            tezhengV{i}=find(lastN(:,i)>=mean(lastN(:,i)));
            tezhengV{i}=sort(tezhengV{i});
            for j=1:length(tezhengV{i})
                tezhengV{i}(j)=mod(tezhengV{i}(j),2^(N-2));
                if  tezhengV{i}(j)==0
                    tezhengV{i}(j)=2^(N-2);
                end
            end
            
            %             arset{1}=tezhengV{1}(1);
            %             a=1;
            for m=1:length(tezhengV)
                arset{m,1}=tezhengV{m}(1);
                a=1;
                for i=2:length(tezhengV{m})
                    if tezhengV{m}(i)-tezhengV{m}(i-1)==1
                        arset{m,a}(end+1)=tezhengV{m}(i);
                    else
                        a=a+1;
                        arset{m,a}=tezhengV{m}(i);
                    end
                end
                
                for i=1:a
                    n=1;
                    for j=1:length(arset{m,i})
                        %                 arset{i}(j)=mod(arset{i}(j),2^(N-2));
                        %                 if  arset{i}(j)==0
                        %                     arset{i}(j)=2^(N-2);
                        %                 end
                        
                        position{m,i}(n)=(1+(arset{m,i}(j)-1)*2)/(2^(N-1));
                        n=n+1;
                    end
                    tz(m,i)=mean(position{m,i});
                end
            end
            
            
            
            %             for i=2:length(tezhengV)
            %                 if tezhengV(i)-tezhengV(i-1)==1
            %                     arset{a}(end+1)=tezhengV(i);
            %                 else
            %                     a=a+1;
            %                     arset{a}=tezhengV(i);
            %                 end
            %             end
            %             for i=1:a
            %                 m=1;
            %                 for j=1:length(arset{i})
            %                     %                 arset{i}(j)=mod(arset{i}(j),2^(N-2));
            %                     %                 if  arset{i}(j)==0
            %                     %                     arset{i}(j)=2^(N-2);
            %                     %                 end
            %
            %                     position{i}(m)=(1+(arset{i}(j)-1)*2)/(2^(N-1));
            %                     m=m+1;
            %                 end
            %                 tz(i)=mean(position{i});
            %             end
        end
        
        %         tezhengV=find(lastN>=mean(lastN));
        %         tezhengV=sort(tezhengV);
        %         for i=1:length(tezhengV)
        %             tezhengV(i)=mod(tezhengV(i),2^(N-2));
        %             if  tezhengV(i)==0
        %                 tezhengV(i)=2^(N-2);
        %             end
        %         end
        % tezhengV  = [1 2 9 10 11 12 13 17 18 20 21 23 24];
        %         arset{1}=tezhengV(1);
        %         a=1;
        %         for i=2:length(tezhengV)
        %             if tezhengV(i)-tezhengV(i-1)==1
        %                 arset{a}(end+1)=tezhengV(i);
        %             else
        %                 a=a+1;
        %                 arset{a}=tezhengV(i);
        %             end
        %         end
        %         for i=1:a
        %             m=1;
        %             for j=1:length(arset{i})
        %                 %                 arset{i}(j)=mod(arset{i}(j),2^(N-2));
        %                 %                 if  arset{i}(j)==0
        %                 %                     arset{i}(j)=2^(N-2);
        %                 %                 end
        %
        %                 position{i}(m)=(1+(arset{i}(j)-1)*2)/(2^(N-1));
        %                 m=m+1;
        %             end
        %             tz(i)=mean(position{i});
        %         end
        % [group,i,j] = findv(n,k);
        % position=(2*j-1)/2^(group-1);
        [DR,DL] = VContinuityInfo1(N);
        
        % CList = [1/8,1/8,0;
        %             0.25,0.25,0;
        %             3/8,3/8,0;
        %             0.5,0.5,0;
        %             5/8,5/8,0;
        %             0.75,0.75,0;
        %             7/8,7/8,0;];
        CList=zeros(2^(N-1)-1,3);
        CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
        CList(:,2)=CList(:,1);
        %       CList(16,:)=[];
        CList(ismember(CList(:,1),tz),:)=[];
        p = sum(CList(:,3)+1);  % 约束方程数目
        SegNum = 2^(N-1); % 分段数
        VNum = 2^N;       % 基函数数目
        knots = linspace(0,1,SegNum+1);% 节点向量
        C = zeros(p,VNum);    % 约束矩阵
        csidx = 0;
        for h = 1 : length(CList)           % 给每个约束建立相应的方程组
            x1 = find(knots==CList(h,1));   % 左节点
            x2 = find(knots==CList(h,2));   % 右节点
            for c = 0 : CList(h,3)          % 从C^0约束至C^CList(h,3)约束，每个约束建立一个方程
                csidx = csidx + 1;
                C(csidx,:) = [DR(:,x1,c+1) - DL(:,x2,c+1)]';
            end
        end
        
        M = [2*A'*A,C';
            C,zeros(p)];
        d = zeros(p,1);
        bx = [2*A'*t;d];
        by = [2*A'*f;d];
        
        X = M\bx;
        Y = M\by;
        
        X = X(1:VNum);
        Y = Y(1:VNum);
        
end

% figure
% bar(abs(Y))
% 找到最大的lamda

if N>2
    %     matrixLamda=ones(length(Y)/2,N);
    %     absY=abs(Y);
    %     matrixLamda(:,1)=kron(absY(1:(k+1)),ones(length(absY)/((k+1)*2),1))/max(absY(1:(k+1)));
    %     for i=2:N
    %         matrixLamda(:,i)=kron(absY(1+(k+1)*2^(i-2):(k+1)*2^(i-1)),ones(length(absY)/((k+1)*2^(i-1)),1))/max(absY(1+(k+1)*2^(i-2):(k+1)*2^(i-1)));
    %     end
    
    %     matrixLamda=ones(length(Y)/2,N);
    %     absY=abs(Y);
    %     matrixLamda(:,1)=kron(absY(1:(k+1)),ones(length(absY)/((k+1)*2),1))/max(absY(1:(k+1)));
    %     for i=2:N
    %         matrixLamda(:,i)=kron(absY(1+(k+1)*2^(i-2):(k+1)*2^(i-1)),ones(length(absY)/((k+1)*2^(i-1)),1))/max(absY(1+(k+1)*2^(i-2):(k+1)*2^(i-1)));
    %     end
    matrixLamda=ones(length(Y)/2,N);
    absY=abs(Y);
    for i=2:N
        for j=1:k+1
            absY((k+1)*2^(i-2)+(j-1)*2^(i-2)+1:(k+1)*2^(i-2)+j*2^(i-2))=absY((k+1)*2^(i-2)+(j-1)*2^(i-2)+1:(k+1)*2^(i-2)+j*2^(i-2))/max(absY((k+1)*2^(i-2)+(j-1)*2^(i-2)+1:(k+1)*2^(i-2)+j*2^(i-2)));
        end
    end
    matrixLamda(:,1)=kron(absY(1:(k+1)),ones(length(absY)/((k+1)*2),1))/max(absY(1:(k+1)));
    for i=2:N
        matrixLamda(:,i)=kron(absY(1+(k+1)*2^(i-2):(k+1)*2^(i-1)),ones(length(absY)/((k+1)*2^(i-1)),1));
    end
    %    0 figure
    %     xlable=[2 N];
    %     ylable=[0.5 4.5];
    %     imagesc(xlable,ylable,matrixLamda(1:length(Y)/2,2:N))
    %     colormap(flipud(bone));
    %
    %     colorbar
    for i=1:k+1
        figure
        xlable=[0 1];
        ylable=[2 N];
        %     colormat
        imagesc(xlable,ylable,matrixLamda(1+length(absY)/(2*(k+1))*(i-1):length(absY)/(2*(k+1))*i,2:N)')
        colormap(flipud(bone));
        %     yticks(2:N)
        set(gca,'YDir','normal');
        colorbar
        caxis([0 1]);
        xlabel('x')
        ylabel('V-系统基函数中的第N组')
        box off
    end
end
%
function [rgx,rgy]=jianruijida(gpoint,k,t)
N=floor(log2(length(gpoint)/(k+1)))+1;
A=LSMatrix_V(k,N,t);
Lambda=A\gpoint;
DC3=Lambda((k+3)*2^(N-2)+1:(k+4)*2^(N-2),:);
[~,locationx] = findpeaks(abs(DC3(:,1)));
[~,locationy] = findpeaks(abs(DC3(:,2)));
xa=(locationx-1)/2^(N-2);
xb=locationx/2^(N-2);
ya=(locationy-1)/2^(N-2);
yb=locationy/2^(N-2);
rgx=[xa,xb];
rgy=[ya,yb];
end

% InfoV1Bas = BaseGene_V1(N);   % 1次V系统基函数信息(非离散采样)
% [VRInfo,NumSeg] = VReconstruction_Polyline([X Y],InfoV1Bas);
% figure,
% plot(t,f,'.','Color',[255 102 102]/255,'MarkerSize',10);hold on
% for j = 1 : NumSeg
%     plot(VRInfo(j,5:6),VRInfo(j,7:8),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
% end
function [phi,K]=qulv(gpoint,i,num)
r=i+1;
l=i-1;
if i==1
    l=num;
elseif i==num
    r=1;
end
v1=gpoint(i,:)-gpoint(l,:);
v2=gpoint(r,:)-gpoint(i,:);
v3=gpoint(r,:)-gpoint(l,:);
phi=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
K=2*sin(phi)/norm(v3);
end
