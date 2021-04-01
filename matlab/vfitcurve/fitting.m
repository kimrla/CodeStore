clear
example=2;
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
        x=0:0.05:10-0.05;        
        y=@(x) 100./exp(abs(x-5))+(x-5).^5/500;
        wide=1;
        for i=1:length(x)
            f(i)=1/wide*integral(y,x(i)-wide/2,x(i)+wide/2);
        end
        xmax=10;
        
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
end

InfoV1Bas = BaseGene_V1(N);   % 1次V系统基函数信息(非离散采样)
[VRInfo,NumSeg] = VReconstruction_Polyline([X Y],InfoV1Bas);
figure,
plot(t,f,'.','Color',[255 102 102]/255,'MarkerSize',10);hold on
for j = 1 : NumSeg
    plot(VRInfo(j,5:6),VRInfo(j,7:8),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
end
xlim([0 1])
% figure
% bar(abs(Y))
% %找到最大的lamda
% n=22;
% [group,i,j] = findv(n,k);
% position=(j)/2^(group-2);
