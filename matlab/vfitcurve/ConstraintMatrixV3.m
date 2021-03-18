function M = ConstraintMatrixV3(N,CList)
%% 3次V-系统的连续性约束矩阵的构造
 % N：V系统层数；CList:节点连续性约束向量m*2(列1：节点值[1..2^(N-1)-1]/列2：连续性值0/1/2/3)
 
k = 3;
SegNum = 2^(N-1);      % 分段数
VNum = (k+1)*SegNum;   % 基函数数目
M = zeros([],VNum);    % 约束矩阵

%% 第2~4项Legendre多项式各阶导数(0~3)信息：(0+)-(1-)
LegInfo = [2*sqrt(3),0,2*sqrt(7);0,-12*sqrt(5),0;0,0,120*sqrt(7);0,0,0];

%% 3次V系统四个生成元函数的各阶导数(0~3)信息：(0+),(1/2+)-(1/2-),(1-)
GenInfo = zeros(4,3,4);
GenInfo(:,:,1) = sqrt(7)*[1,0,1;-18,0,18;132,0,132;-384,768,384];
GenInfo(:,:,2) = sqrt(5)*[1,0,-1;-24,0,-24;228,384,-228;-840,0,-840];
GenInfo(:,:,3) = sqrt(3)*[1,0,1;-28,80,28;312,0,312;-1344,2688,1344];
GenInfo(:,:,4) = [1,8,-1;-30,0,-30;360,960,-360;-1680,0,-1680];

%%
PNum = size(CList,1);  % 约束点数目
p = 2^(N-1);           % 约束点位置分母
csidx = 0;             % 约束方程指标
for h = 1 : PNum
    q = CList(h,1);    % 约束点h位置分子
    
    if q == 0          % 端点处连续性约束（对Legendre多项式有约束）
       ch = CList(h,2);% 端点处对应的连续性0/1/2/3 
       for c = 0 : ch
           csidx = csidx + 1;
           M(csidx,2:4) = LegInfo(c+1,:);   % Legendre基约束
           for s = 5 : 8                    % 生成元基约束
               M(csidx,s) = GenInfo(c+1,1,s-4) - GenInfo(c+1,3,s-4);
           end
           for n = 3 : N                    % 后续各层基函数约束
               j1 = 1;
               j2 = 2^(n-2);
               for i = 1 : 4
                   s1 = (k+i)*2^(n-2)+j1;
                   s2 = (k+i)*2^(n-2)+j2;
                   M(csidx,s1) = sqrt(2^(n-2))*GenInfo(c+1,1,i);
                   M(csidx,s2) = -sqrt(2^(n-2))*GenInfo(c+1,3,i);
               end
           end
       end
    else               % 内部节点处连续性约束
    g = gcd(p,q);
    ph = p / g;        % 约分后约束点h位置分母
    qh = q / g;
    n = log2(ph)+1;    % 约束点h对应的层数和序数
    j = (qh+1)/2;
    
    % 根据n和ch确定需约束的生成元基函数类；n层之后每层每类基函数有2项受约束
    ch = CList(h,2);    % 约束点h对应的连续性0/1/2/3
    for c = 0 : ch
        csidx = csidx + 1;
        
        for i = (4-ch) : 4           % 第n层受约束基函数
            s = (k+i)*2^(n-2) + j;   % 受约束基函数的序数
            M(csidx,s) = 2^(c*(n-2))*sqrt(2^(n-2))*GenInfo(c+1,2,i);
        end
        
        for L = n+1 : N      % 第n+1,n+2,...,N层受约束基函数
            j1 = (L-n)*(2*j-1);    % 该层约束基函数左
            %j2 = j1 + 1;          % 该层约束基函数右
            for i = 1 : 4
                s1 = (k+i)*2^(L-2)+j1;
                s2 = s1 + 1;
                M(csidx,s1) = -2^(c*(L-2))*sqrt(2^(L-2)) * GenInfo(c+1,3,i);
                M(csidx,s2) = 2^(c*(L-2))*sqrt(2^(L-2)) * GenInfo(c+1,1,i);
            end
        end
    end
    end
end

    
end