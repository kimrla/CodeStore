
function [PP_all,T_all,TP_all] = PIA_CurApr_TUpdate(P,k,pro,NumItr)
%% 曲线拟合的PIA方法(更新控制点后重新参数化)
%% InPut: P-型值点列N*2
%%        k-B样条次数
%%        pro-型值延拓类型：‘1’-端点插值；‘2’-周期延拓
%%        NumItr:迭代次数

%% 数据准备
n = k + 1;                     % 阶数
N = size(P,1);                 % 型值点数目（若为封闭点列，首末端点应相同）

PP_all = zeros(N+2,NumItr+1);  % 迭代各步(延拓)控制点列
T_all = zeros(N+6,NumItr+1);   % 迭代各步对应的基函数节点向量
TP_all = zeros(N+2,NumItr+1);
[PP,T,TP] = ParaP(P,k,pro);
PP_all(:,1) = complex(PP(:,1),PP(:,2)); 
T_all(:,1) = T;
TP_all(:,1) = TP;

%% 迭代
for m = 1 : NumItr
    
    P_0 = [real(PP_all(:,m)),imag(PP_all(:,m))];  % 更新前延拓控制点列
    T_0 = T_all(:,m);                             % 更新前基函数节点矢量
    TP_0 = TP_all(:,m);
    
    P_1 = P_0 + PP - bspline_deboor(n,T_0,P_0,TP_0);
    
    [P_2,T_all(:,m+1),TP_all(:,m+1)] = ParaP(P_1(2:end-1,:),k,pro);
    
    PP_all(:,m+1) = complex(P_2(:,1),P_2(:,2));

end

%% 
 function  [PP,T,TP] = ParaP(P,k,pro)
 %% 计算给定型值点列的延拓（PP）、基函数节点（T）、型值参数（TP）
    %% 对型值做弦长参数化
    Chord = vecnorm(diff(P,1,1),2,2);                                % 型值点列各段弦长(共N-1段)
    NormChord = [0;cumsum(Chord)/sum(Chord)];                        % 各型值的归一化弦长参数
    %% 基函数节点矢量
    if k == 2                            
       T0 = (NormChord(1:end-1) + NormChord(2:end))/2;
    elseif k == 3
       T0 = NormChord; 
    end
    
    switch pro
    case 1  % 开放型曲线
        PP = [P(1,:);P;P(end,:)];    % 型值延拓
%         PP_all(:,1) = complex(PP(:,1),PP(:,2));
        T = [0;0;0;T0;1;1;1];        % 基函数节点矢量
        TP = [0;NormChord;1];        % 延拓后型值参数
      
    case 2  % 周期型曲线
        PP = [P(end-1,:);P;P(2,:)];  % 型值延拓
%         PP_all(:,1) = complex(PP(:,1),PP(:,2));
        if k==2
           T = [T0(end-2:end)-1;T0;T0(1:3)+1];  
        elseif k==3
           T = [T0(end-3:end-1)-1;T0;T0(2:4)+1];   
        end 
        TP = [NormChord(end-1);NormChord;NormChord(2)];  % 延拓后型值参数
    end
    
 end

end