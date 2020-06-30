function [N,R,P] = kongzhidingdian(M,n,p,t,ui,d)
%KONGZHIDINGDIAN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
R=zeros(1,2);
for i=2:M   
    for j = 2 : n
        N(i, j) = Njp(j, p , t(i), ui);%Nij�洢Njp(ti) i=1~M+1,j=0~n
    end
    R(i,:)=d(i,:)-Njp(1,p,t(i),ui)*d(1,:)-Njp(n+1,p,t(i),ui)*d(M,:);% R(i)=d(i)-N0pti*d(1)-Nnpti*d(M);
end
P=pinv(transpose(N)*N)*transpose(N)*R;%NT*N*P=NT*R;������ƶ���Pi i=1~n+1 
P(1,:)=d(1,:);%������P0=d0 Pn=dM  -> P1=d1 Pn+1=dM+1 
P(n+1,:)=d(M+1,:);
end

