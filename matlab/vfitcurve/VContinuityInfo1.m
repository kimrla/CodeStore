function [DR,DL] = VContinuityInfo1(N)
%% 1次V-系统前N组[(k+1)2^(N -1)项]基函数在节点处的连续性信息

k = 1;
VNum = (k+1)*2^(N-1);             % 基函数数目
knots = linspace(0,1,2^(N-1)+1);  % 节点向量

DR = zeros(VNum,2^(N-1)+1,k+1);   % 基函数在各节点处的0~k阶右导数
DL = DR;                          % 基函数在各节点处的0~k阶左导数

%% Legendre多项式(1次项)
l2 = sqrt(3)*[1-2*knots;-2*ones(size(knots))]; % Legendre多项式(1次项)的0阶与1阶导数值
for i = 1 : 2
    DR(2,1:end-1,i) = l2(i,1:end-1);   
    DL(2,2:end,i) = l2(i,2:end);
end

%% 生成元（3~4）
%% 右导数
x1 = knots(find(knots<1/2));
x2 = knots(find(knots>=1/2 & knots<1));
% 0阶导数值
DR(3:4,1:end-1,1) = [sqrt(3)*[-4 * x1 + 1, 4 * x2 - 3];
                             [-6 * x1 + 1,-6 * x2 + 5]];
% 1阶导数值
DR(3:4,1:end-1,2) = [sqrt(3)*[repmat(-4,[1,size(x1,2)]),repmat(4,[1,size(x2,2)])];
                    [repmat(-6,[1,size(x1,2)]),repmat(-6,[1,size(x2,2)])]]; 

%% 左导数                
x1 = knots(find(knots>0 & knots<=1/2));
x2 = knots(find(knots>1/2));
% 0阶导数值
DL(3:4,2:end,1) = [sqrt(3)*[-4 * x1 + 1, 4 * x2 - 3];
                           [-6 * x1 + 1,-6 * x2 + 5]];
% 1阶导数值                 
DL(3:4,2:end,2) = [sqrt(3)*[repmat(-4,[1,size(x1,2)]),repmat(4,[1,size(x2,2)])];
                           [repmat(-6,[1,size(x1,2)]),repmat(-6,[1,size(x2,2)])]];                

DR_Gen = DR(3:4,1:end-1,:);
DL_Gen = DL(3:4,2:end,:);  

%% 后续各层
for n = 3 : N
    NumJ = 2^(n-2);                                 % 第n层每类基函数数目
    KontSeg = [reshape([1:2^(N-1)]',[],NumJ)]'; 
%     KontSeg_R = [reshape([1:2^(N-1)]',[],NumJ)]';
%     KontSeg_L = [reshape([2:2^(N-1)+1]',[],NumJ)]';
    DR_n = sqrt(NumJ)*DR_Gen(:,1:NumJ:end,:);
    DL_n = sqrt(NumJ)*DL_Gen(:,fliplr(end:-NumJ:1),:);
    for j = 1 : NumJ
        vidx = [0;1]*NumJ+j+2^(n-1);
        DR(vidx,KontSeg(j,:),:) = DR_n;
        DL(vidx,KontSeg(j,:)+1,:) = DL_n;
    end
 
    % 每求一阶导数, 系数乘以2^(n-2) 
    DR(2^(n-1)+1:2^n,:,2) = 2^(n-2)*DR(2^(n-1)+1:2^n,:,2);
    DL(2^(n-1)+1:2^n,:,2) = 2^(n-2)*DL(2^(n-1)+1:2^n,:,2);
      
end
                
                
% %% Legendre多项式(2~4)的0~3阶导数值
% function [g2,g3,g4] = LegDerivatives(x)
%     g2 = sqrt(3)*[1-2*x;-2*ones(size(x));0*ones(size(x));0*ones(size(x))];
%     g3 = sqrt(5)*[6 * x.^2 - 6 * x + 1;12 * x - 6;12*ones(size(x));0*ones(size(x))];
%     g4 = sqrt(7)*[-20 * x.^3 + 30 * x.^2 - 12 * x + 1;-60 * x.^2 + 60 * x - 12;-120 * x + 60;-120*ones(size(x))];
%     %g = [g2 g3 g4];
% end

% %% 生成元
% function g = GenDerivatives(x)
% 
% if x>=0 && x<1/2
%    g1 = sqrt(7)*[-64 * x.^3 + 66 * x.^2 - 18 * x + 1;-192 * x.^2 + 132 * x - 18;-384 * x + 132;-384]; 
%    g2 = sqrt(5)*[-140 * x.^3 + 114 * x.^2 - 24 * x + 1;-420 * x.^2 + 228 * x - 24;-840 * x + 228;-840];
%    g3 = sqrt(3)*[-224 * x.^3 + 156 * x.^2 - 28 * x + 1;-672 * x.^2 + 312 * x - 28;-1344 * x + 312;-1344];
%    g4 = [-280 * x.^3 + 180 * x.^2 - 30 * x + 1;-840 * x.^2 + 360 * x - 30;-1680 * x + 360;-1680];
%    g = [g1 g2 g3 g4];
% elseif x>1/2 && x<=1
%     g1 = sqrt(7)*[64 * x.^3 - 126 * x.^2 + 78 * x - 15;192 * x.^2 - 252 * x + 78;384 * x - 252;384];
%     g2 = sqrt(5)*[-140 * x.^3 + 306 * x.^2 - 216 * x + 49;-420 * x.^2 + 612 * x - 216;-840 * x + 612;-840];
%     g3 = sqrt(3)*[224 * x.^3 - 516 * x.^2 + 388 * x - 95;672 * x.^2 - 1032 * x + 388;1344 * x - 1032;1344];
%     g4 = [-280 * x.^3 + 660 * x.^2 - 510 * x + 129;-840 * x.^2 + 1320 * x - 510;-1680*x+1320;-1680];
%     g = [g1 g2 g3 g4];
% elseif x == 1/2
%     g1 = sqrt(7)*[0;0;0;768];
%     g2 = sqrt(5)*[0;0;384;0];
%     g3 = sqrt(3)*[0;80;0;2688];
%     g4 = [8;0;960;0];
%     g = [g1 g2 g3 g4];
% else 
%     error('ERROR!!输入的x不在[0,1]内！');
% end
% 
% end


end