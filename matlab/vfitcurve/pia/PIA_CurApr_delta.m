function [PP,T] = PIA_CurApr_delta(P,k,pro)
%% ������ϵ�PIA����(��������ֵdelta)
 % ÿ�ε���ֻ������������ֵ����ֵ��

%% ����׼��
n = k + 1;                                                       % ����
N = size(P,1);                                                   % ��ֵ����Ŀ����Ϊ��յ��У���ĩ�˵�Ӧ��ͬ��

nItrMax = 10;   % ������������

Chord = vecnorm(diff(P,1,1),2,2);                                % ��ֵ���и����ҳ�(��N-1��)
NormChord = [0;cumsum(Chord)/sum(Chord)];                        % ����ֵ�Ĺ�һ���ҳ�����

if k == 2                            
   T0 = (NormChord(1:end-1) + NormChord(2:end))/2;
elseif k == 3
   T0 = NormChord; 
end

switch pro
    case 1  % ����������
        PP = [P(1,:);P;P(end,:)];    % ��ֵ����
        T = [0;0;0;T0;1;1;1];        % �������ڵ�ʸ��
        TP = [0;NormChord;1];        % ���غ���ֵ����
      
    case 2  % ����������
        PP = [P(end-1,:);P;P(2,:)];  % ��ֵ����
        if k==2
           T = [T0(end-2:end)-1;T0;T0(1:3)+1];  
        elseif k==3
           T = [T0(end-3:end-1)-1;T0;T0(2:4)+1];   
        end 
        TP = [NormChord(end-1);NormChord;NormChord(2)];  % ���غ���ֵ����
end

%% ����
m = 0;
P_0 = P;
while m < nItrMax
    
end


end