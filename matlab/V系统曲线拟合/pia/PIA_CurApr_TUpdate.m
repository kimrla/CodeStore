
function [PP_all,T_all,TP_all] = PIA_CurApr_TUpdate(P,k,pro,NumItr)
%% ������ϵ�PIA����(���¿��Ƶ�����²�����)
%% InPut: P-��ֵ����N*2
%%        k-B��������
%%        pro-��ֵ�������ͣ���1��-�˵��ֵ����2��-��������
%%        NumItr:��������

%% ����׼��
n = k + 1;                     % ����
N = size(P,1);                 % ��ֵ����Ŀ����Ϊ��յ��У���ĩ�˵�Ӧ��ͬ��

PP_all = zeros(N+2,NumItr+1);  % ��������(����)���Ƶ���
T_all = zeros(N+6,NumItr+1);   % ����������Ӧ�Ļ������ڵ�����
TP_all = zeros(N+2,NumItr+1);
[PP,T,TP] = ParaP(P,k,pro);
PP_all(:,1) = complex(PP(:,1),PP(:,2)); 
T_all(:,1) = T;
TP_all(:,1) = TP;

%% ����
for m = 1 : NumItr
    
    P_0 = [real(PP_all(:,m)),imag(PP_all(:,m))];  % ����ǰ���ؿ��Ƶ���
    T_0 = T_all(:,m);                             % ����ǰ�������ڵ�ʸ��
    TP_0 = TP_all(:,m);
    
    P_1 = P_0 + PP - bspline_deboor(n,T_0,P_0,TP_0);
    
    [P_2,T_all(:,m+1),TP_all(:,m+1)] = ParaP(P_1(2:end-1,:),k,pro);
    
    PP_all(:,m+1) = complex(P_2(:,1),P_2(:,2));

end

%% 
 function  [PP,T,TP] = ParaP(P,k,pro)
 %% ���������ֵ���е����أ�PP�����������ڵ㣨T������ֵ������TP��
    %% ����ֵ���ҳ�������
    Chord = vecnorm(diff(P,1,1),2,2);                                % ��ֵ���и����ҳ�(��N-1��)
    NormChord = [0;cumsum(Chord)/sum(Chord)];                        % ����ֵ�Ĺ�һ���ҳ�����
    %% �������ڵ�ʸ��
    if k == 2                            
       T0 = (NormChord(1:end-1) + NormChord(2:end))/2;
    elseif k == 3
       T0 = NormChord; 
    end
    
    switch pro
    case 1  % ����������
        PP = [P(1,:);P;P(end,:)];    % ��ֵ����
%         PP_all(:,1) = complex(PP(:,1),PP(:,2));
        T = [0;0;0;T0;1;1;1];        % �������ڵ�ʸ��
        TP = [0;NormChord;1];        % ���غ���ֵ����
      
    case 2  % ����������
        PP = [P(end-1,:);P;P(2,:)];  % ��ֵ����
%         PP_all(:,1) = complex(PP(:,1),PP(:,2));
        if k==2
           T = [T0(end-2:end)-1;T0;T0(1:3)+1];  
        elseif k==3
           T = [T0(end-3:end-1)-1;T0;T0(2:4)+1];   
        end 
        TP = [NormChord(end-1);NormChord;NormChord(2)];  % ���غ���ֵ����
    end
    
 end

end