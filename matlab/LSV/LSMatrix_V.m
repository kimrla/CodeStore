function A = LSMatrix_V(k,N,T)
%% k次V-系统最小二乘曲线拟合的系数矩阵
%% N ：所用V-系统层数；T：数据点对应参数向量(列)

%% Check
if min(T)<0 || max(T)>1
    error('the parameters is beyond the range [0,1]');
end

%%
M = length(T);   
NumV = (k+1)*2^(N-1);
A = zeros(M,NumV);
A(:,1:k+1) = LegendreValue(k,T);       % Legendre多项式
A(:,k+2:2*(k+1)) = GenerateValue(k,T); % 生成元
for n = 3 : N 
    NumSeg = 2^(n-2);
    knots = linspace(0,1,NumSeg+1);
    for j = 1 : NumSeg
        t0 = knots(j);
        t1 = knots(j+1);
        Tidx = find(T>=t0 & T<=t1);    % 子区间T索引          
        Tj = NumSeg*(T(Tidx) - t0);    % 变换至[0,1]区间
        
        Gj = sqrt(NumSeg)*GenerateValue(k,Tj);   % 本层第j组基函数
        NumPreV = (k+1)*NumSeg;
        Vidx = [(NumPreV + j):NumSeg:(2*NumPreV)];
        A(Tidx,Vidx) = Gj;
         
        if (~isempty(find(T==t0))) && (j~=1)
            A(find(T==t0),Vidx) = A(find(T==t0),Vidx)/2;
        end
        
        if (~isempty(find(T==t1))) && (j~=NumSeg)
            A(find(T==t1),Vidx) = A(find(T==t1),Vidx)/2;
        end
    end
end

    function L = LegendreValue(k,x)
    % Legendre多项式
       L = ones(length(x),4);
       L(:,2) = sqrt(3)*(1-2*x);
       L(:,3) = sqrt(5)*(6 * x.^2 - 6 * x + 1);
       L(:,4) = sqrt(7)*(-20 * x.^3 + 30 * x.^2 - 12 * x + 1);
       
       L = L(:,1:k+1);
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

end