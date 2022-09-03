function [PP,T] = PIA_CurApr_delta(P,k,pro)
%% 曲线拟合的PIA方法(拟合误差阈值delta)
 % 每次迭代只更新误差大于阈值的型值点

%% 数据准备
n = k + 1;                                                       % 阶数
N = size(P,1);                                                   % 型值点数目（若为封闭点列，首末端点应相同）

nItrMax = 10;   % 迭代次数上限

Chord = vecnorm(diff(P,1,1),2,2);                                % 型值点列各段弦长(共N-1段)
NormChord = [0;cumsum(Chord)/sum(Chord)];                        % 各型值的归一化弦长参数

if k == 2                            
   T0 = (NormChord(1:end-1) + NormChord(2:end))/2;
elseif k == 3
   T0 = NormChord; 
end

switch pro
    case 1  % 开放型曲线
        PP = [P(1,:);P;P(end,:)];    % 型值延拓
        T = [0;0;0;T0;1;1;1];        % 基函数节点矢量
        TP = [0;NormChord;1];        % 延拓后型值参数
      
    case 2  % 周期型曲线
        PP = [P(end-1,:);P;P(2,:)];  % 型值延拓
        if k==2
           T = [T0(end-2:end)-1;T0;T0(1:3)+1];  
        elseif k==3
           T = [T0(end-3:end-1)-1;T0;T0(2:4)+1];   
        end 
        TP = [NormChord(end-1);NormChord;NormChord(2)];  % 延拓后型值参数
end

%% 迭代
m = 0;
P_0 = P;
while m < nItrMax
    
end


end