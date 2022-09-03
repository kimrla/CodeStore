function [P,NumSeg] = VReconstruction_Polyline(VS,V1)
%% 折线的V1系统重构
 % VS：V谱(N*2);V1:1次V系统信息矩阵(非离散采样)
 
 NumV = length(VS);      % V谱数目
 NumLev = log2(NumV);
 if(NumLev ~= fix(NumLev))
    error('错误！输入V谱数目应为2^n!');
 end
 
 NumSeg = 2^(NumLev-1);        % 重构折线段数
 NumBas = 2 * NumSeg;
 
 P = zeros(NumSeg,8);          % 重构折线信息，每段保存斜率和截距  左右两端函数值[x(left),x(right),y(left),y(right)]
 
 t = linspace(0,1,NumSeg+1);
%  sq3 = sqrt(3);
%  P(:,2) = VS(1,1);
%  P(:,4) = VS(1,2);
%  
%  P(:,1:2) = repmat(VS(2,1)*sq3*[-2,1],[NumSeg,1]);
%  P(:,3:4) = repmat(VS(2,2)*sq3*[-2,1],[NumSeg,1]);
 
 %% 计算斜率和截距
  for i = 1 : NumBas
     ti = V1(i,1:3);                                        % 基函数i的支集节点(t1,t2,t3)
     td = [find(t==ti(1)),find(t==ti(2)),find(t==ti(3))];   % 支集节点索引
     
     vi = [V1(i,4:5);V1(i,6:7)];                               % 基函数i左右两段
     NumHlf = td(2) - td(1);
     for s = 1 : 2
         P(td(s):td(s+1)-1,1:2) = P(td(s):td(s+1)-1,1:2) + repmat(VS(i,1)*vi(s,:),[NumHlf,1]);  % x
         P(td(s):td(s+1)-1,3:4) = P(td(s):td(s+1)-1,3:4) + repmat(VS(i,2)*vi(s,:),[NumHlf,1]);  % y
     end    
  end
  
  %% 计算节点函数值
  for j = 1 : NumSeg
      tj = [t(j:j+1);1 1];
      P(j,5:6) = P(j,1:2)*tj;
      P(j,7:8) = P(j,3:4)*tj;
  end
  
  %% 画图
%   figure,
%   for j = 1 : NumSeg
%       plot(P(j,5:6),P(j,7:8),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
%   end
%   axis equal
 
end