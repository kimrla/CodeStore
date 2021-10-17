
function [PP_all,T] = PIA_CurApr(P,k,pro,NumItr)
%% 曲线拟合的PIA方法
%% InPut: P-型值点列N*2
%%        k-B样条次数
%%        pro-型值延拓类型：‘1’-端点插值；‘2’-周期延拓
%%        NumItr:迭代次数

%% 数据准备
n = k + 1;                                                       % 阶数
N = size(P,1);                                                   % 型值点数目（若为封闭点列，首末端点应相同）

Chord = vecnorm(diff(P,1,1),2,2);                                % 型值点列各段弦长(共N-1段)
NormChord = [0;cumsum(Chord)/sum(Chord)];                        % 各型值的归一化弦长参数

% exd_left = ceil((k-1)/2);   % 型值向左延拓数目
% exd_rigt = floor((k-1)/2);  % 型值向右延拓数目

PP_all = zeros(N+2,NumItr+1);        % 迭代各步(延拓)控制点列

%% 基函数节点矢量
if k == 2                            
   T0 = (NormChord(1:end-1) + NormChord(2:end))/2;
elseif k == 3
   T0 = NormChord; 
end

switch pro
    case 1  % 开放型曲线
        PP = [P(1,:);P;P(end,:)];    % 型值延拓
        PP_all(:,1) = complex(PP(:,1),PP(:,2));
        T = [0;0;0;T0;1;1;1];        % 基函数节点矢量
        TP = [0;NormChord;1];        % 延拓后型值参数
      
    case 2  % 周期型曲线
        PP = [P(end-1,:);P;P(2,:)];  % 型值延拓
        PP_all(:,1) = complex(PP(:,1),PP(:,2));
        if k==2
           T = [T0(end-2:end)-1;T0;T0(1:3)+1];  
        elseif k==3
           T = [T0(end-3:end-1)-1;T0;T0(2:4)+1];   
        end 
        TP = [NormChord(end-1);NormChord;NormChord(2)];  % 延拓后型值参数
end

%% 迭代
%  PP_1 = PP;
for m = 1 : NumItr
    PP_0 = [real(PP_all(:,m)),imag(PP_all(:,m))];                      % 更新前延拓控制点列
    PP_1 = PP_0 + PP - bspline_deboor(n,T,PP_0,TP);
    PP_all(:,m+1) = complex(PP_1(:,1),PP_1(:,2));
end
% 
% switch pro
%     case 1  % 开放型曲线
%         T = [zeros(k,1);NormChord;ones(k,1)];                                 % 节点矢量
%         PP = [repmat(P(1,:),[exd_left,1]);P;repmat(P(end,:),[exd_rigt,1])];   % 型值延拓
%         PP_all(:,1) = complex(PP(:,1),PP(:,2));
%         for m = 1 : NumItr
%             PP_0 = [real(PP_all(:,m)),imag(PP_all(:,m))];                      % 更新前延拓控制点列
%             PP_1 = PP_0(exd_left+1:end-exd_rigt,:) + P - bspline_deboor(n,T,PP_0,NormChord);
%             PP_1 = [repmat(PP_1(1,:),[exd_left,1]);PP_1;repmat(PP_1(end,:),[exd_rigt,1])];   % 型值延拓
%             PP_all(:,m+1) = complex(PP_1(:,1),PP_1(:,2));
%         end
%         
%     case 2  % 封闭型曲线
%         T = [NormChord(end-k:end-1)-1;NormChord;NormChord(2:k+1)+1];       % 节点矢量
%         PP = [P(end-exd_left:end-1,:);P;P(2:1+exd_rigt,:)];                % 型值延拓
%         PP_all(:,1) = complex(PP(:,1),PP(:,2));
%         for m = 1 : NumItr
%             PP_0 = [real(PP_all(:,m)),imag(PP_all(:,m))];                      % 更新前延拓控制点列
%             PP_1 = PP_0(exd_left+1:end-exd_rigt,:) + P - bspline_deboor(n,T,PP_0,NormChord);
%             PP_1 = [PP_1(end-exd_left:end-1,:);PP_1;PP_1(2:1+exd_rigt,:)];   % 型值延拓
%             PP_all(:,m+1) = complex(PP_1(:,1),PP_1(:,2));
%         end
% end

% T_PP = T(k-exd_left+1:end-k+exd_rigt);   % 延拓型值点对应参数值
% U = linspace(0,1,2001);                  % B样条曲线采样点对应参数值

%% 型值点(控制点)更新
% for m = 1 : NumItr
%     PP_0 = [real(PP_all(:,m)),imag(PP_all(:,m))];  % 更新前控制点列
%     PP_1 = PP_0 + PP - bspline_deboor(n,T,PP_0,T_PP);
%     PP_all(:,m+1) = complex(PP_1(:,1),PP_1(:,2));
% end

% C0 = bspline_deboor(n,T,[real(PP_all(:,1)),imag(PP_all(:,1))],U);
% C1 = bspline_deboor(n,T,[real(PP_all(:,end)),imag(PP_all(:,end))],U);
% figure,
% plot(P(:,1),P(:,2),'.','MarkerSize',10);hold on
% plot(P(:,1),P(:,2));hold on
% plot(C0(:,1),C0(:,2));hold on
% plot(C1(:,1),C1(:,2));
% %plot(C_T(:,1),C_T(:,2),'.','MarkerSize',15);hold on
% axis equal

end