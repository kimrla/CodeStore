
%% 
load V3_4096  % O

%% Test1
x = [1/8192 : 1/4096 : 1]';
N = 4;
SegNum = 2^(N-1); % ·Ö¶ÎÊý
VNum = 2^(N+1);


CList = [4,3];
M = ConstraintMatrixV3(N,CList);
rM = rank(M);

Z = null(M);

beta = rand(size(Z,2),1);

f = O(1:VNum,:)' * (Z * beta);
xx = reshape(x,[],SegNum);
ff = reshape(f,[],SegNum);
figure,
for s = 1 : SegNum
    plot(xx(:,s),ff(:,s),'LineWidth',2);hold on
end