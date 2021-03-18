
%% ��1:��Բ��������1��V����Լ����С���˱ƽ�

clear all
close all

%% Test1
t = [1/64:1/32:1]';
x = cos(2*pi*t);
y = sin(2*pi*t);
% figure,
% plot(x,y,'.','MarkerSize',10);
% axis equal

k = 1;
N = 4;

%% Լ������
[DR,DL] = VContinuityInfo1(N);
% CList = [0,1,0;
%          0.5,0.5,0;
%          0.25,0.25,0;
%          0.75,0.75,0];
CList = [0.5,0.5,0;
         0.25,0.25,0;
         0.75,0.75,0;
         1/8,1/8,0;
         3/8,3/8,0;
         5/8,5/8,0;
         7/8,7/8,0;
         0,1,0];
p = sum(CList(:,3)+1);  % Լ��������Ŀ
SegNum = 2^(N-1); % �ֶ���
VNum = 2^N;       % ��������Ŀ
knots = linspace(0,1,SegNum+1);% �ڵ�����
C = zeros(p,VNum);    % Լ������
csidx = 0; 
for h = 1 : length(CList)           % ��ÿ��Լ��������Ӧ�ķ�����
    x1 = find(knots==CList(h,1));   % ��ڵ�
    x2 = find(knots==CList(h,2));   % �ҽڵ�
    for c = 0 : CList(h,3)          % ��C^0Լ����C^CList(h,3)Լ����ÿ��Լ������һ������
        csidx = csidx + 1;
        C(csidx,:) = [DR(:,x1,c+1) - DL(:,x2,c+1)]';
    end
end  

%% ϵ������
A = LSMatrix_V(k,N,t);

M = [2*A'*A,C';
     C,zeros(p)];
d = zeros(p,1);
bx = [2*A'*x;d]; 
by = [2*A'*y;d]; 

X = M\bx;
Y = M\by;

X = X(1:VNum);
Y = Y(1:VNum);

%% Show
InfoV1Bas = BaseGene_V1(N);   % 1��Vϵͳ��������Ϣ(����ɢ����)
[VRInfo,NumSeg] = VReconstruction_Polyline([X Y],InfoV1Bas);
figure,
plot(x,y,'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
for j = 1 : NumSeg
    plot(VRInfo(j,5:6),VRInfo(j,7:8),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
end
axis equal
axis off