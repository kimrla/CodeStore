
function [PP_all,T] = PIA_CurApr(P,k,pro,NumItr)
%% ������ϵ�PIA����
%% InPut: P-��ֵ����N*2
%%        k-B��������
%%        pro-��ֵ�������ͣ���1��-�˵��ֵ����2��-��������
%%        NumItr:��������

%% ����׼��
n = k + 1;                                                       % ����
N = size(P,1);                                                   % ��ֵ����Ŀ����Ϊ��յ��У���ĩ�˵�Ӧ��ͬ��

Chord = vecnorm(diff(P,1,1),2,2);                                % ��ֵ���и����ҳ�(��N-1��)
NormChord = [0;cumsum(Chord)/sum(Chord)];                        % ����ֵ�Ĺ�һ���ҳ�����

% exd_left = ceil((k-1)/2);   % ��ֵ����������Ŀ
% exd_rigt = floor((k-1)/2);  % ��ֵ����������Ŀ

PP_all = zeros(N+2,NumItr+1);        % ��������(����)���Ƶ���

%% �������ڵ�ʸ��
if k == 2                            
   T0 = (NormChord(1:end-1) + NormChord(2:end))/2;
elseif k == 3
   T0 = NormChord; 
end

switch pro
    case 1  % ����������
        PP = [P(1,:);P;P(end,:)];    % ��ֵ����
        PP_all(:,1) = complex(PP(:,1),PP(:,2));
        T = [0;0;0;T0;1;1;1];        % �������ڵ�ʸ��
        TP = [0;NormChord;1];        % ���غ���ֵ����
      
    case 2  % ����������
        PP = [P(end-1,:);P;P(2,:)];  % ��ֵ����
        PP_all(:,1) = complex(PP(:,1),PP(:,2));
        if k==2
           T = [T0(end-2:end)-1;T0;T0(1:3)+1];  
        elseif k==3
           T = [T0(end-3:end-1)-1;T0;T0(2:4)+1];   
        end 
        TP = [NormChord(end-1);NormChord;NormChord(2)];  % ���غ���ֵ����
end

%% ����
%  PP_1 = PP;
for m = 1 : NumItr
    PP_0 = [real(PP_all(:,m)),imag(PP_all(:,m))];                      % ����ǰ���ؿ��Ƶ���
    PP_1 = PP_0 + PP - bspline_deboor(n,T,PP_0,TP);
    PP_all(:,m+1) = complex(PP_1(:,1),PP_1(:,2));
end
% 
% switch pro
%     case 1  % ����������
%         T = [zeros(k,1);NormChord;ones(k,1)];                                 % �ڵ�ʸ��
%         PP = [repmat(P(1,:),[exd_left,1]);P;repmat(P(end,:),[exd_rigt,1])];   % ��ֵ����
%         PP_all(:,1) = complex(PP(:,1),PP(:,2));
%         for m = 1 : NumItr
%             PP_0 = [real(PP_all(:,m)),imag(PP_all(:,m))];                      % ����ǰ���ؿ��Ƶ���
%             PP_1 = PP_0(exd_left+1:end-exd_rigt,:) + P - bspline_deboor(n,T,PP_0,NormChord);
%             PP_1 = [repmat(PP_1(1,:),[exd_left,1]);PP_1;repmat(PP_1(end,:),[exd_rigt,1])];   % ��ֵ����
%             PP_all(:,m+1) = complex(PP_1(:,1),PP_1(:,2));
%         end
%         
%     case 2  % ���������
%         T = [NormChord(end-k:end-1)-1;NormChord;NormChord(2:k+1)+1];       % �ڵ�ʸ��
%         PP = [P(end-exd_left:end-1,:);P;P(2:1+exd_rigt,:)];                % ��ֵ����
%         PP_all(:,1) = complex(PP(:,1),PP(:,2));
%         for m = 1 : NumItr
%             PP_0 = [real(PP_all(:,m)),imag(PP_all(:,m))];                      % ����ǰ���ؿ��Ƶ���
%             PP_1 = PP_0(exd_left+1:end-exd_rigt,:) + P - bspline_deboor(n,T,PP_0,NormChord);
%             PP_1 = [PP_1(end-exd_left:end-1,:);PP_1;PP_1(2:1+exd_rigt,:)];   % ��ֵ����
%             PP_all(:,m+1) = complex(PP_1(:,1),PP_1(:,2));
%         end
% end

% T_PP = T(k-exd_left+1:end-k+exd_rigt);   % ������ֵ���Ӧ����ֵ
% U = linspace(0,1,2001);                  % B�������߲������Ӧ����ֵ

%% ��ֵ��(���Ƶ�)����
% for m = 1 : NumItr
%     PP_0 = [real(PP_all(:,m)),imag(PP_all(:,m))];  % ����ǰ���Ƶ���
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