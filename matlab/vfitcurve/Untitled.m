clear,clc
% 
% xiaoshu2dec(0.1,5);
% 
% function y=xiaoshu2dec(x,N)
% y=zeros(1,N);
% for i=1:N    
%     y(i)=floor(x*2);
%     x=x*2-y(i);
% end
% end
N = 6;
x = 1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
k = 1;

% A = LSMatrix_V(k,N,x);
a=[1,2,3;4,5,6;0,8,9];
C=orth(a);
B=Gram_Schmidt_Orthogonalization(a);
% B=Schmidt_orthogonalization(A);
E=B*B';
F=C'*C;
% function b=Schmidt_orthogonalization(a)
% [m,n] = size(a);
% if(m<n)
%     error('行小于列，无法计算，请转置后重新输入');
%     return
% end
% b=zeros(m,n);
% %正交化
% b(:,1)=a(:,1);
% for i=2:n
%     for j=1:i-1
%         b(:,i)=b(:,i)-dot(a(:,i),b(:,j))/dot(b(:,j),b(:,j))*b(:,j);
%     end
%     b(:,i)=b(:,i)+a(:,i);
% end
% end
function T = Gram_Schmidt_Orthogonalization(p_tr)

    % 一列为一个向量

    [row,col]= size(p_tr);

    T = zeros(row,col);

    

    T(:,1)=p_tr(:,1);

    for i = 2 : col

        for j = 1: i-1

            p_tr(:,i)= p_tr(:,i) - ((T(:,j)' * p_tr(:,i))/(T(:,j)' * T(:,j))) * T(:,j);

        end

        T(:,i)=p_tr(:,i);

    end

    % 向量单位化

    for i = 1: col

        length=norm(T(:,i));

        for j = 1: row

            T(j,i)= T(j,i)/ length;

        end

    end

  end