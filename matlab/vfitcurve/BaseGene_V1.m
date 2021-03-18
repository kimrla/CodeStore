function V = BaseGene_V1(NumLev)
%% 1��Vϵͳ������(����ɢ����)
 % NumLev������(Legendre����ʽΪ��1��)
 % V: ÿ��������Ϊ2��1�ζ���ʽ������֧���ڵ�(x1,x2,x3),ÿ��1�ζ���ʽ��б�ʡ��ؾࡢ���˺���ֵ��11������
 
 NumBase = 2^NumLev;     % ��������Ŀ
 V = zeros(NumBase,11);  % ��������Ϣ[x1,x2,x3,k_left,b_left,k_right,b_right,v_left(x1),v_left(x2),v_right(x2),v_right(x3)]
                         % ��������Ϣ[x1,x2,x3,k_left,b_left,k_right,b_right,v_left(x1),v_left(x2),v_right(x2),v_right(x3)]
 sq3 = sqrt(3);
 
 V(1,:) = [0,0.5,1,0,1,0,1,1,1,1,1];
 V(2,:) = [0,0.5,1,-2*sq3,sq3,-2*sq3,sq3,sq3,0,0,-sq3];       % Legendre����ʽ����2��
 
%  V(3,:) = [0,0.5,1,-4*sqrt(3),sqrt(3),4*sqrt(3),-3*sqrt(3)];  % ����Ԫ
%  V(4,:) = [0,0.5,1,-6,1,-6,5];
 
 for L = 2 : NumLev
     NumPre = 2^(L-1);  % ǰL-1���������Ŀ
     NumThs = NumPre;   % �����������Ŀ
     Num_i = NumThs/2;  % ����ÿ���������Ŀ
     n2 = 2^L;
     c1 = sqrt(n2/2);
     c2 = sqrt(n2/4); 
     for j = 1 : Num_i
         x1 = (j-1)/Num_i;
         x3 = j/Num_i;
         x2 = (x1+x3)/2;
         V(NumPre+j,1:3) = [x1,x2,x3];
         V(NumPre+Num_i+j,1:3) = [x1,x2,x3];
         
         V(NumPre+j,4:7) = sq3*c2*[-n2,4*j-3,n2,-4*j+1];
         V(NumPre+Num_i+j,4:7) = c2*[-1.5*n2,6*j-5,-1.5*n2,6*j-1];
         
         V(NumPre+j,8:11) = sq3*c2*[1,-1,-1,1];
         V(NumPre+Num_i+j,8:11) = [c2,-c1,c1,-c2];
     end
 end

end