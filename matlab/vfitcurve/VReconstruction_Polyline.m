function [P,NumSeg] = VReconstruction_Polyline(VS,V1)
%% ���ߵ�V1ϵͳ�ع�
 % VS��V��(N*2);V1:1��Vϵͳ��Ϣ����(����ɢ����)
 
 NumV = length(VS);      % V����Ŀ
 NumLev = log2(NumV);
 if(NumLev ~= fix(NumLev))
    error('��������V����ĿӦΪ2^n!');
 end
 
 NumSeg = 2^(NumLev-1);        % �ع����߶���
 NumBas = 2 * NumSeg;
 
 P = zeros(NumSeg,8);          % �ع�������Ϣ��ÿ�α���б�ʺͽؾ�  �������˺���ֵ[x(left),x(right),y(left),y(right)]
 
 t = linspace(0,1,NumSeg+1);
%  sq3 = sqrt(3);
%  P(:,2) = VS(1,1);
%  P(:,4) = VS(1,2);
%  
%  P(:,1:2) = repmat(VS(2,1)*sq3*[-2,1],[NumSeg,1]);
%  P(:,3:4) = repmat(VS(2,2)*sq3*[-2,1],[NumSeg,1]);
 
 %% ����б�ʺͽؾ�
  for i = 1 : NumBas
     ti = V1(i,1:3);                                        % ������i��֧���ڵ�(t1,t2,t3)
     td = [find(t==ti(1)),find(t==ti(2)),find(t==ti(3))];   % ֧���ڵ�����
     
     vi = [V1(i,4:5);V1(i,6:7)];                               % ������i��������
     NumHlf = td(2) - td(1);
     for s = 1 : 2
         P(td(s):td(s+1)-1,1:2) = P(td(s):td(s+1)-1,1:2) + repmat(VS(i,1)*vi(s,:),[NumHlf,1]);  % x
         P(td(s):td(s+1)-1,3:4) = P(td(s):td(s+1)-1,3:4) + repmat(VS(i,2)*vi(s,:),[NumHlf,1]);  % y
     end    
  end
  
  %% ����ڵ㺯��ֵ
  for j = 1 : NumSeg
      tj = [t(j:j+1);1 1];
      P(j,5:6) = P(j,1:2)*tj;
      P(j,7:8) = P(j,3:4)*tj;
  end
  
  %% ��ͼ
%   figure,
%   for j = 1 : NumSeg
%       plot(P(j,5:6),P(j,7:8),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
%   end
%   axis equal
 
end