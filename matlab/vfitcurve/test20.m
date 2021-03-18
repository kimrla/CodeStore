
%% ��������һ�����ߣ��������(C^0����)������ڵ�ΪC^2����

clear all

load V3_4096  % O

%% Test1
x = [1/8192 : 1/4096 : 1]';
N = 3;
SegNum = 2^(N-1);                % �ֶ���
VNum = 2^(N+1);                  % ��������Ŀ

[DR,DL] = VContinuityInfo3(N);   % �������ڽڵ㴦���Ҹ��׵���
knots = linspace(0,1,2^(N-1)+1); % �ڵ������������׽ڵ�x=0,ĩ�ڵ�x=1���Լ��м�ڵ㣩

%% �ڵ㴦������Լ������
CList = [1/4,1/4,0;   % �Ҳ༫�ޣ���༫�ޣ�������
         1/2,1/2,2;
         3/4,3/4,2];
     
M = zeros([],VNum);    % Լ������
csidx = 0; 
for h = 1 : size(CList,1)
    x1 = find(knots==CList(h,1));   
    x2 = find(knots==CList(h,2));
    for c = 0 : CList(h,3)
        csidx = csidx + 1;
        M(csidx,:) = [DR(:,x1,c+1) - DL(:,x2,c+1)]';
    end 
end     

Z = null(M);            % Z ���пռ�ΪԼ���ӿռ���������

beta = 10 * rand(size(Z,2),1);  

f = O(1:VNum,:)' * (Z * beta);
xx = reshape(x,[],SegNum);
ff = reshape(f,[],SegNum);
figure,
for s = 1 : SegNum
    plot(xx(:,s),ff(:,s),'LineWidth',2);hold on
end


