function Lambda = LSCurFit_VPIA(P,k,N,t,CList,j)
%% ����V-ϵͳ��(����)Լ����С�������߱ƽ�
 % P:�������(NumP*NumDim)��k:������N��������������CList:�ڵ�������Լ��
 % ˼·��1)�г�������Ԫ���������нڵ㴦��2^(N-1)����0,1,...,k���ҵ���[0,1)������(0,1]
 %       2)CList��ÿ��Ϊ[��ڵ�,�ҽڵ�,������]��ʽ���ֽ�������л���������ڵ㴦�ҵ������ҽڵ㴦����
 %         ����ڵ�Ϊ��,
%  %% ���в�����(�ҳ�������)
 [NumP,NumDim] = size(P);
%  csl = [0;cumsum(vecnorm(diff(P),2,2).^0.5)];       % �ۼӳ���
%  LChord = csl(end)+norm(P(1,:)-P(end,:))^0.5;
%  t = csl / LChord;
%  

 %% ϵ������
 A = LSMatrix_V(k,N,t);
 A=A(:,size(A,2)/2+j:2^(N-2):end);
% A=A(utidx,size(A,2)/2+1:end);
 %% Լ������
 NumCns = sum(CList(:,3)+1);   % Լ��������Ŀ
 NumBas = (k+1)*2^(N-1);       % ��������Ŀ
 C = zeros(NumCns,NumBas);     % Լ������
 cidx = 1; 
 for h = 1 : size(CList,1)
     C(cidx:cidx+CList(h,3),:) = derivativediff(k,N,CList(h,1),CList(h,2),CList(h,3));  % ����0~CList(h,3)�׵���֮��
     cidx = cidx + CList(h,3) + 1;
 end
 C=C(:,size(C,2)/2+j:2^(N-2):end); 
% C=C(:,size(C,2)/2+1:end); 
 %% KKT
%  NumBas=NumBas/2;
 M = [2*A'*A,C';
      C,zeros(NumCns)];               % KKT����ϵ������
 b = [2*A'*P;zeros(NumCns,NumDim)] ;  % 
%  
% M=[A;C];
% b=[P;zeros(NumCns,NumDim)];

%  Lambda0 = M \ b;
Lambda0=pinv(M)*b;
 Lambda = Lambda0(1:size(A,2),:);
%  Lambda = Lambda0(1:NumBas,:);
end

%%
function C = derivativediff(k,N,x1,x2,m)
%% ����k��V-ϵͳǰN��������ڽڵ�x1�ҵ�����ڵ�x2����֮�0��,1��,...,m��
NumBas = (k+1)*2^(N-1);       % ��������Ŀ
C = zeros(m+1,NumBas);

%% Legendre����ʽ
 % ���г�ȫ��0~3��Legendre����ʽ��x1(x2)����0~3�׵���
L1 = [[1;0;0;0],...
       sqrt(3)*[1-2*x1;-2;0;0],...
       sqrt(5)*[6*x1^2-6*x1+1;12*x1-6;12;0],...
       sqrt(7)*[-20*x1^3+30*x1^2-12*x1+1;-60*x1^2+60*x1-12;-120*x1+60;-120]];
L2 = [[1;0;0;0],...
       sqrt(3)*[1-2*x2;-2;0;0],...
       sqrt(5)*[6*x2^2-6*x2+1;12*x2-6;12;0],...
       sqrt(7)*[-20*x2^3+30*x2^2-12*x2+1;-60*x2^2+60*x2-12;-120*x2+60;-120]];   

C(:,1:k+1) = L1(1:m+1,1:k+1) - L2(1:m+1,1:k+1);

%% ����Ԫ
 % �г�k+1��k������Ԫ��x1��0~k���ҵ�����x2��0~m������
C(:,k+2:2*k+2) = GenDiff(k,x1,x2,m);

%% ����������
for n = 3 : N
    numb = 2^(n-2);              % ����ÿ���������Ŀ
    X = linspace(0,1,numb+1);
    for i = 1 : numb
        if x1>=X(i) && x1<X(i+1)
           j1 = i;
           break
        end
    end
    for i = 1 : numb
        if x2>X(i) && x2<=X(i+1)
           j2 = i;
           break
        end
    end
    % ��x1,x2����������Ԫ��Ӧ�ڵ�
    xx1 = numb*x1 - j1 + 1;
    xx2 = numb*x2 - j2 + 1;
    
    scl = sqrt(numb)*(numb.^[0:m]');  % m�׵���ϵ��2^[m*(n-2)]
    scls = repmat(scl,1,k+1);
    
    NumPre = (k+1)*2^(n-2);         % 1~(n-1)���������Ŀ
    J = numb*[0:k];
    if j1 == j2 
       C(:,NumPre+J+j1) = scls .* GenDiff(k,xx1,xx2,m); 
    else
       C(:,NumPre+J+j1) = scls .* GenDiff(k,xx1,0,m);
       C(:,NumPre+J+j2) = scls .* GenDiff(k,1,xx2,m);
    end
    
end
 
%% 
function D = GenDiff(k,x1,x2,m)
% k��V-ϵͳ����Ԫ�����ڽڵ�x1��0~m���ҵ��� �� �ڽڵ�x2��0~m������֮��
% D : (m+1)*(k+1)
if x1 > x2
    error('ERROR!');
end
if m > k
    error('ERROR!');
end
 %%%%% x1
x = x1;
if x>=0 && x<1/2
   switch k
       case 1
           G = [sqrt(3)*[-4*x+1;-4],...
                        [-6*x+1;-6]];
       case 2
           G = [sqrt(5)*[16*x^2-10*x+1;32*x-10;32],...
                sqrt(3)*[30*x^2-14*x+1;60*x-14;60],...
                        [40*x^2-16*x+1;80*x-16;80]];
       case 3
           G = [sqrt(7)*[-64*x^3+66*x^2-18*x+1;-192*x^2+132*x-18;-384*x+132;-384],...
                sqrt(5)*[-140*x^3+114*x^2-24*x+1;-420*x^2+228*x-24;-840*x+228;-840],...
                sqrt(3)*[-224 * x^3 + 156 * x^2 - 28 * x + 1;-672 * x^2 + 312 * x - 28;-1344 * x + 312;-1344],...
                        [-280 * x^3 + 180 * x^2 - 30 * x + 1;-840 * x^2 + 360 * x - 30;-1680 * x + 360;-1680]];
           
   end
elseif x>=1/2 && x<1
   switch k
       case 1
           G = [sqrt(3)*[4*x-3;4],...
                        [-6*x+5;-6]];
       case 2
           G = [sqrt(5)*[-16*x^2+22*x-7;-32*x+22;-32],...
                sqrt(3)*[30*x^2-46*x+17;60*x-46;60],...
                        [-40*x^2+64*x-25;-80*x+64;-80]];
       case 3
           G = [sqrt(7)*[64 * x^3 - 126 * x^2 + 78 * x - 15;192 * x^2 - 252 * x + 78;384 * x - 252;384],...
                sqrt(5)*[-140 * x^3 + 306 * x^2 - 216 * x + 49;-420 * x^2 + 612 * x - 216;-840 * x + 612;-840],...
                sqrt(3)*[224 * x^3 - 516 * x^2 + 388 * x - 95;672 * x^2 - 1032 * x + 388;1344 * x - 1032;1344],...
                        [-280 * x^3 + 660 * x^2 - 510 * x + 129;-840 * x^2 + 1320 * x - 510;-1680*x+1320;-1680]];     
   end 
else
    %error('ERROR!the knot x1 is not in [0,1)');
    G = zeros(k+1);
end

 %%%%%%% x2
x = x2;
if x>0 && x<=1/2
   switch k
       case 1
           H = [sqrt(3)*[-4*x+1;-4],...
                        [-6*x+1;-6]];
       case 2
           H = [sqrt(5)*[16*x^2-10*x+1;32*x-10;32],...
                sqrt(3)*[30*x^2-14*x+1;60*x-14;60],...
                        [40*x^2-16*x+1;80*x-16;80]];
       case 3
           H = [sqrt(7)*[-64*x^3+66*x^2-18*x+1;-192*x^2+132*x-18;-384*x+132;-384],...
                sqrt(5)*[-140*x^3+114*x^2-24*x+1;-420*x^2+228*x-24;-840*x+228;-840],...
                sqrt(3)*[-224 * x^3 + 156 * x^2 - 28 * x + 1;-672 * x^2 + 312 * x - 28;-1344 * x + 312;-1344],...
                        [-280 * x^3 + 180 * x^2 - 30 * x + 1;-840 * x^2 + 360 * x - 30;-1680 * x + 360;-1680]];
           
   end
elseif x>1/2 && x<=1
   switch k
       case 1
           H = [sqrt(3)*[4*x-3;4],...
                        [-6*x+5;-6]];
       case 2
           H = [sqrt(5)*[-16*x^2+22*x-7;-32*x+22;-32],...
                sqrt(3)*[30*x^2-46*x+17;60*x-46;60],...
                        [-40*x^2+64*x-25;-80*x+64;-80]];
       case 3
           H = [sqrt(7)*[64 * x^3 - 126 * x^2 + 78 * x - 15;192 * x^2 - 252 * x + 78;384 * x - 252;384],...
                sqrt(5)*[-140 * x^3 + 306 * x^2 - 216 * x + 49;-420 * x^2 + 612 * x - 216;-840 * x + 612;-840],...
                sqrt(3)*[224 * x^3 - 516 * x^2 + 388 * x - 95;672 * x^2 - 1032 * x + 388;1344 * x - 1032;1344],...
                        [-280 * x^3 + 660 * x^2 - 510 * x + 129;-840 * x^2 + 1320 * x - 510;-1680*x+1320;-1680]];     
   end 
else
    %error('ERROR!the knot x2 is not in (0,1]');
    H = zeros(k+1);
end

 % x1 - x2
D = G(1:m+1,:) - H(1:m+1,:);  
end


end
