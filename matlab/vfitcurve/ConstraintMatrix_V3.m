function M = ConstraintMatrix_V3(N,CList)
%% 3��V-ϵͳ��������Լ������Ĺ���
 % N��Vϵͳ������CList:�ڵ�������Լ������m*3(��1/2���ڵ�ֵ[1..2^(N-1)-1]/��3��������ֵ0/1/2/3)
 
k = 3;
SegNum = 2^(N-1);      % �ֶ���
VNum = (k+1)*SegNum;   % ��������Ŀ
M = zeros([],VNum);    % Լ������

%% ��2~4��Legendre����ʽ���׵���(0~3)��Ϣ��(0+)-(1-)
LegInfo = [2*sqrt(3),0,2*sqrt(7);0,-12*sqrt(5),0;0,0,120*sqrt(7);0,0,0];

%% 3��Vϵͳ�ĸ�����Ԫ�����ĸ��׵���(0~3)��Ϣ��(0+),(1/2+)-(1/2-),(1-)
GenInfo = zeros(4,3,4);
GenInfo(:,:,1) = sqrt(7)*[1,0,1;-18,0,18;132,0,132;-384,768,384];
GenInfo(:,:,2) = sqrt(5)*[1,0,-1;-24,0,-24;228,384,-228;-840,0,-840];
GenInfo(:,:,3) = sqrt(3)*[1,0,1;-28,80,28;312,0,312;-1344,2688,1344];
GenInfo(:,:,4) = [1,8,-1;-30,0,-30;360,960,-360;-1680,0,-1680];

%%
PNum = size(CList,1);  % Լ������Ŀ 
p = 2^(N-1);           % Լ����λ�÷�ĸ
csidx = 0;             % Լ������ָ��
for h = 1 : PNum
    %% q1~=q2
    ch = CList(h,3);
    cc = [0:ch]';
    Mm = zeros(ch+1,VNum,2);
    for m = 1 : 2
        q = CList(h,m);
        x = q/p;
        
        %% Legendre��
        v = LegDerivatives(x);  % 4*3
        v = v(1:ch+1,:);
        Mm(:,2:4,m) = v;
        
       %% Ѱ��x���ڲ�n��j
        g = gcd(p,q);
        ph = p / g;        % Լ�ֺ�Լ����hλ�÷�ĸ
        qh = q / g;
        n = log2(ph)+1;    % Լ����h��Ӧ�Ĳ���������
        j = (qh+1)/2;
        
        for L = 2 : (n-1)
            % x��L���е�jֵ
            knots = linspace(0,1,2^(L-2)+1);
            jL = find((knots - x)>0,1)-1;
            
        end
        
    end
    
    q1 = CList(h,1);
    q2 = CList(h,2);
    x1 = q1 / p;
    x2 = q2 / p;
    ch = CList(h,3);    % Լ����h��Ӧ��������0/1/2/3
    Mh = zeros(ch+1,VNum);

    % n = 1��Legendre����ʽ��
       v1 = LegDerivatives(x1);  % 4*3
       v1 = v1(1:ch+1,:);
       v2 = LegDerivatives(x2);
       v2 = v2(1:ch+1,:);
       
       Mh(:,2:4) = (v1 - v2); 
    end
    
    for n = 2 : N
        J = 2^(n-2);
        
    end
    
    for c = 0 : ch
        
    end
    
    
    %% ������������
    if q1 == q2
       q = q1;
       g = gcd(p,q);
       ph = p / g;        % Լ�ֺ�Լ����hλ�÷�ĸ
       qh = q / g;
       n = log2(ph)+1;    % Լ����h��Ӧ�Ĳ���������
       j = (qh+1)/2;
       for c = 0 : ch
           csidx = csidx + 1;
           % ��n����Լ��������
           for i = (4-ch) : 4            
               s = (k+i)*2^(n-2) + j;    % ��Լ��������������
               M(csidx,s) = 2^(c*(n-2))*sqrt(2^(n-2))*GenInfo(c+1,2,i);
           end
           % ��n+1,n+2,...,N����Լ��������
           for L = n+1 : N      
              j1 = (L-n)*(2*j-1);        % �ò�Լ����������
               for i = 1 : 4
                s1 = (k+i)*2^(L-2)+j1;   % Լ�������������
                s2 = s1 + 1;             % Լ���һ���������
                M(csidx,s1) = -2^(c*(L-2))*sqrt(2^(L-2)) * GenInfo(c+1,3,i);
                M(csidx,s2) = 2^(c*(L-2))*sqrt(2^(L-2)) * GenInfo(c+1,1,i);
               end
           end
       end
    %% ˫������Լ��  
    else
        for m = 1 : 2          % ������(m=1ʱ�ҵ�����m=2ʱ����)
            if (q2 - q1) < 2
                error('ERROR!!');
            end
            q = CList(h,m);
            g = gcd(p,q);
           ph = p / g;         % Լ�ֺ�Լ����hλ�÷�ĸ
           qh = q / g;
            n = log2(ph)+1;    % Լ����h��Ӧ�Ĳ���������
            j = (qh+1)/2;
            for c = 0 : ch
                csidx = csidx + 1;
                
            end
        end
    end
    
end

%% 3��V-ϵͳ����Ԫ0~3�׵���
function g = GenDerivatives(x)

if x>=0 && x<1/2
   g1 = sqrt(7)*[-64 * x^3 + 66 * x^2 - 18 * x + 1;-192 * x^2 + 132 * x - 18;-384 * x + 132;-384]; 
   g2 = sqrt(5)*[-140 * x^3 + 114 * x^2 - 24 * x + 1;-420 * x^2 + 228 * x - 24;-840 * x + 228;-840];
   g3 = sqrt(3)*[-224 * x^3 + 156 * x^2 - 28 * x + 1;-672 * x^2 + 312 * x - 28;-1344 * x + 312;-1344];
   g4 = [-280 * x^3 + 180 * x^2 - 30 * x + 1;-840 * x^2 + 360 * x - 30;-1680 * x + 360;-1680];
   g = [g1 g2 g3 g4];
elseif x>1/2 && x<=1
    g1 = sqrt(7)*[64 * x^3 - 126 * x^2 + 78 * x - 15;192 * x^2 - 252 * x + 78;384 * x - 252;384];
    g2 = sqrt(5)*[-140 * x^3 + 306 * x^2 - 216 * x + 49;-420 * x^2 + 612 * x - 216;-840 * x + 612;-840];
    g3 = sqrt(3)*[224 * x^3 - 516 * x^2 + 388 * x - 95;672 * x^2 - 1032 * x + 388;1344 * x - 1032;1344];
    g4 = [-280 * x^3 + 660 * x^2 - 510 * x + 129;-840 * x^2 + 1320 * x - 510;-1680*x+1320;-1680];
    g = [g1 g2 g3 g4];
elseif x == 1/2
    g1 = sqrt(7)*[0;0;0;768];
    g2 = sqrt(5)*[0;0;384;0];
    g3 = sqrt(3)*[0;80;0;2688];
    g4 = [8;0;960;0];
    g = [g1 g2 g3 g4];
else 
    error('ERROR!!�����x����[0,1]�ڣ�');
end

end

%% 3��V-ϵͳLegendre����ʽ0~3�׵���
function g = LegDerivatives(x)
    g2 = sqrt(3)*[1-2*x;-2;0;0];
    g3 = sqrt(5)*[6 * x^2 - 6 * x + 1;12 * x - 6;12;0];
    g4 = sqrt(7)*[-20 * x^3 + 30 * x^2 - 12 * x + 1;-60 * x^2 + 60 * x - 12;-120 * x + 60;-120];
    g = [g2 g3 g4];
end
    
end