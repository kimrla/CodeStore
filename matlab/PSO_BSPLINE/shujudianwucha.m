function [epsilon,e] = shujudianwucha(M,N,P,d)
%SHUJUDIANWUCHA �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
epsilon=zeros(1,2);
e=0;
for i=1:M+1
        %����i=di-sti=di-��Njpti*Pj
        s(i,:)=N(i,:)*P;
        epsilon(i,:)=d(i,:)-N(i,:)*P; 
        e=e+norm(epsilon(i,:))^2;%�������
end  
end

