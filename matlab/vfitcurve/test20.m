
%% 例：构造一条曲线，带尖锐点(C^0连续)，其余节点为C^2连续

clear all

load V3_4096  % O

%% Test1
x = [1/8192 : 1/4096 : 1]';
N = 3;
SegNum = 2^(N-1);                % 分段数
VNum = 2^(N+1);                  % 基函数数目

[DR,DL] = VContinuityInfo3(N);   % 基函数在节点处左右各阶导数
knots = linspace(0,1,2^(N-1)+1); % 节点向量（包括首节点x=0,末节点x=1，以及中间节点）

%% 节点处连续性约束条件
CList = [1/4,1/4,0;   % 右侧极限，左侧极限，连续性
         1/2,1/2,2;
         3/4,3/4,2];
     
M = zeros([],VNum);    % 约束矩阵
csidx = 0; 
for h = 1 : size(CList,1)
    x1 = find(knots==CList(h,1));   
    x2 = find(knots==CList(h,2));
    for c = 0 : CList(h,3)
        csidx = csidx + 1;
        M(csidx,:) = [DR(:,x1,c+1) - DL(:,x2,c+1)]';
    end 
end     

Z = null(M);            % Z 的列空间为约束子空间中正交基

beta = 10 * rand(size(Z,2),1);  

f = O(1:VNum,:)' * (Z * beta);
xx = reshape(x,[],SegNum);
ff = reshape(f,[],SegNum);
figure,
for s = 1 : SegNum
    plot(xx(:,s),ff(:,s),'LineWidth',2);hold on
end


