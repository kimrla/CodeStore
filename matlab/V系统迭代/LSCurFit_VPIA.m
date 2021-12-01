function Lambda = LSCurFit_VPIA(P,k,N,t,CList,j)
%% 基于V-系统的(线性)约束最小二乘曲线逼近
 % P:有序点列(NumP*NumDim)；k:次数；N：基函数组数；CList:节点连续性约束
 % 思路：1)列出各生成元函数在所有节点处（2^(N-1)）的0,1,...,k阶右导数[0,1)和左导数(0,1]
 %       2)CList中每行为[左节点,右节点,连续性]形式，分解算出所有基函数在左节点处右导数与右节点处左导数
 %         以左节点为例,
%  %% 点列参数化(弦长参数化)
 [NumP,NumDim] = size(P);
%  csl = [0;cumsum(vecnorm(diff(P),2,2).^0.5)];       % 累加长度
%  LChord = csl(end)+norm(P(1,:)-P(end,:))^0.5;
%  t = csl / LChord;
%  

 %% 系数矩阵
 A = LSMatrix_V(k,N,t);
 A=A(:,size(A,2)/2+j:2^(N-2):end);
% A=A(utidx,size(A,2)/2+1:end);
 %% 约束方程
 NumCns = sum(CList(:,3)+1);   % 约束方程数目
 NumBas = (k+1)*2^(N-1);       % 基函数数目
 C = zeros(NumCns,NumBas);     % 约束矩阵
 cidx = 1; 
 for h = 1 : size(CList,1)
     C(cidx:cidx+CList(h,3),:) = derivativediff(k,N,CList(h,1),CList(h,2),CList(h,3));  % 计算0~CList(h,3)阶导数之差
     cidx = cidx + CList(h,3) + 1;
 end
 C=C(:,size(C,2)/2+j:2^(N-2):end); 
% C=C(:,size(C,2)/2+1:end); 
 %% KKT
%  NumBas=NumBas/2;
 M = [2*A'*A,C';
      C,zeros(NumCns)];               % KKT方程系数矩阵
 b = [2*A'*P;zeros(NumCns,NumDim)] ;  % 
%  
% M=[A;C];
% b=[P;zeros(NumCns,NumDim)];

%  Lambda0 = M \ b;
Lambda0=pinv(M)*b;
 Lambda = Lambda0(1:size(A,2),:);
%  Lambda = Lambda0(1:NumBas,:);
end

%%
function C = derivativediff(k,N,x1,x2,m)
%% 计算k次V-系统前N组基函数在节点x1右导数与节点x2左导数之差，0阶,1阶,...,m阶
NumBas = (k+1)*2^(N-1);       % 基函数数目
C = zeros(m+1,NumBas);

%% Legendre多项式
 % 先列出全部0~3次Legendre多项式在x1(x2)处的0~3阶导数
L1 = [[1;0;0;0],...
       sqrt(3)*[1-2*x1;-2;0;0],...
       sqrt(5)*[6*x1^2-6*x1+1;12*x1-6;12;0],...
       sqrt(7)*[-20*x1^3+30*x1^2-12*x1+1;-60*x1^2+60*x1-12;-120*x1+60;-120]];
L2 = [[1;0;0;0],...
       sqrt(3)*[1-2*x2;-2;0;0],...
       sqrt(5)*[6*x2^2-6*x2+1;12*x2-6;12;0],...
       sqrt(7)*[-20*x2^3+30*x2^2-12*x2+1;-60*x2^2+60*x2-12;-120*x2+60;-120]];   

C(:,1:k+1) = L1(1:m+1,1:k+1) - L2(1:m+1,1:k+1);

%% 生成元
 % 列出k+1个k次生成元在x1处0~k阶右导数与x2处0~m阶左导数
C(:,k+2:2*k+2) = GenDiff(k,x1,x2,m);

%% 后续基函数
for n = 3 : N
    numb = 2^(n-2);              % 本层每类基函数数目
    X = linspace(0,1,numb+1);
    for i = 1 : numb
        if x1>=X(i) && x1<X(i+1)
           j1 = i;
           break
        end
    end
    for i = 1 : numb
        if x2>X(i) && x2<=X(i+1)
           j2 = i;
           break
        end
    end
    % 将x1,x2换算至生成元对应节点
    xx1 = numb*x1 - j1 + 1;
    xx2 = numb*x2 - j2 + 1;
    
    scl = sqrt(numb)*(numb.^[0:m]');  % m阶导数系数2^[m*(n-2)]
    scls = repmat(scl,1,k+1);
    
    NumPre = (k+1)*2^(n-2);         % 1~(n-1)层基函数数目
    J = numb*[0:k];
    if j1 == j2 
       C(:,NumPre+J+j1) = scls .* GenDiff(k,xx1,xx2,m); 
    else
       C(:,NumPre+J+j1) = scls .* GenDiff(k,xx1,0,m);
       C(:,NumPre+J+j2) = scls .* GenDiff(k,1,xx2,m);
    end
    
end
 
%% 
function D = GenDiff(k,x1,x2,m)
% k次V-系统生成元函数在节点x1处0~m阶右导数 与 在节点x2处0~m阶左导数之差
% D : (m+1)*(k+1)
if x1 > x2
    error('ERROR!');
end
if m > k
    error('ERROR!');
end
 %%%%% x1
x = x1;
if x>=0 && x<1/2
   switch k
       case 1
           G = [sqrt(3)*[-4*x+1;-4],...
                        [-6*x+1;-6]];
       case 2
           G = [sqrt(5)*[16*x^2-10*x+1;32*x-10;32],...
                sqrt(3)*[30*x^2-14*x+1;60*x-14;60],...
                        [40*x^2-16*x+1;80*x-16;80]];
       case 3
           G = [sqrt(7)*[-64*x^3+66*x^2-18*x+1;-192*x^2+132*x-18;-384*x+132;-384],...
                sqrt(5)*[-140*x^3+114*x^2-24*x+1;-420*x^2+228*x-24;-840*x+228;-840],...
                sqrt(3)*[-224 * x^3 + 156 * x^2 - 28 * x + 1;-672 * x^2 + 312 * x - 28;-1344 * x + 312;-1344],...
                        [-280 * x^3 + 180 * x^2 - 30 * x + 1;-840 * x^2 + 360 * x - 30;-1680 * x + 360;-1680]];
           
   end
elseif x>=1/2 && x<1
   switch k
       case 1
           G = [sqrt(3)*[4*x-3;4],...
                        [-6*x+5;-6]];
       case 2
           G = [sqrt(5)*[-16*x^2+22*x-7;-32*x+22;-32],...
                sqrt(3)*[30*x^2-46*x+17;60*x-46;60],...
                        [-40*x^2+64*x-25;-80*x+64;-80]];
       case 3
           G = [sqrt(7)*[64 * x^3 - 126 * x^2 + 78 * x - 15;192 * x^2 - 252 * x + 78;384 * x - 252;384],...
                sqrt(5)*[-140 * x^3 + 306 * x^2 - 216 * x + 49;-420 * x^2 + 612 * x - 216;-840 * x + 612;-840],...
                sqrt(3)*[224 * x^3 - 516 * x^2 + 388 * x - 95;672 * x^2 - 1032 * x + 388;1344 * x - 1032;1344],...
                        [-280 * x^3 + 660 * x^2 - 510 * x + 129;-840 * x^2 + 1320 * x - 510;-1680*x+1320;-1680]];     
   end 
else
    %error('ERROR!the knot x1 is not in [0,1)');
    G = zeros(k+1);
end

 %%%%%%% x2
x = x2;
if x>0 && x<=1/2
   switch k
       case 1
           H = [sqrt(3)*[-4*x+1;-4],...
                        [-6*x+1;-6]];
       case 2
           H = [sqrt(5)*[16*x^2-10*x+1;32*x-10;32],...
                sqrt(3)*[30*x^2-14*x+1;60*x-14;60],...
                        [40*x^2-16*x+1;80*x-16;80]];
       case 3
           H = [sqrt(7)*[-64*x^3+66*x^2-18*x+1;-192*x^2+132*x-18;-384*x+132;-384],...
                sqrt(5)*[-140*x^3+114*x^2-24*x+1;-420*x^2+228*x-24;-840*x+228;-840],...
                sqrt(3)*[-224 * x^3 + 156 * x^2 - 28 * x + 1;-672 * x^2 + 312 * x - 28;-1344 * x + 312;-1344],...
                        [-280 * x^3 + 180 * x^2 - 30 * x + 1;-840 * x^2 + 360 * x - 30;-1680 * x + 360;-1680]];
           
   end
elseif x>1/2 && x<=1
   switch k
       case 1
           H = [sqrt(3)*[4*x-3;4],...
                        [-6*x+5;-6]];
       case 2
           H = [sqrt(5)*[-16*x^2+22*x-7;-32*x+22;-32],...
                sqrt(3)*[30*x^2-46*x+17;60*x-46;60],...
                        [-40*x^2+64*x-25;-80*x+64;-80]];
       case 3
           H = [sqrt(7)*[64 * x^3 - 126 * x^2 + 78 * x - 15;192 * x^2 - 252 * x + 78;384 * x - 252;384],...
                sqrt(5)*[-140 * x^3 + 306 * x^2 - 216 * x + 49;-420 * x^2 + 612 * x - 216;-840 * x + 612;-840],...
                sqrt(3)*[224 * x^3 - 516 * x^2 + 388 * x - 95;672 * x^2 - 1032 * x + 388;1344 * x - 1032;1344],...
                        [-280 * x^3 + 660 * x^2 - 510 * x + 129;-840 * x^2 + 1320 * x - 510;-1680*x+1320;-1680]];     
   end 
else
    %error('ERROR!the knot x2 is not in (0,1]');
    H = zeros(k+1);
end

 % x1 - x2
D = G(1:m+1,:) - H(1:m+1,:);  
end


end
