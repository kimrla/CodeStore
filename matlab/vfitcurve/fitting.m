clear
example=1;
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
        x=0:0.02:10-0.02;
        f=100./exp(abs(x-5))+(x-5).^5/500;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        xmax=10;
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
% % end
k = 1;
N = 7;
A = LSMatrix_V(k,N,t);

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
        X = A\t;
        Y=A\f;
        % V = LSMatrix_V(k,N,linspace(0,1,501)');
        % F = V * X;
        % figure,
        % plot(t,f,'.r');hold on
        % plot(linspace(0,1,501)',F,'b');
        
        
        
    case 2
        % 约束方程
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

figure
bar(abs(Y))
% 找到最大的lamda

matrixLamda=ones(length(Y)/2,N);
matrixLamda(:,1)=kron(Y(1:2),ones(32,1));
for i=2:N
matrixLamda(:,i)=kron(Y(1+(k+1)*2^(i-2):(k+1)*2^(i-1)),ones(length(Y)/2^i,1));
end
matrixLamda=abs(matrixLamda);
imagesc(matrixLamda(:,2:N))
colorbar

InfoV1Bas = BaseGene_V1(N);   % 1次V系统基函数信息(非离散采样)
[VRInfo,NumSeg] = VReconstruction_Polyline([X Y],InfoV1Bas);
figure,
plot(t,f,'.','Color',[255 102 102]/255,'MarkerSize',10);hold on
for j = 1 : NumSeg
    plot(VRInfo(j,5:6),VRInfo(j,7:8),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
end
xlim([0 1])
