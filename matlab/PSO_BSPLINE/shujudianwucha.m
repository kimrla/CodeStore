function [epsilon,e] = shujudianwucha(M,N,P,d)
%SHUJUDIANWUCHA 此处显示有关此函数的摘要
%   此处显示详细说明
epsilon=zeros(1,2);
e=0;
for i=1:M+1
        %误差εi=di-sti=di-ΣNjpti*Pj
        s(i,:)=N(i,:)*P;
        epsilon(i,:)=d(i,:)-N(i,:)*P; 
        e=e+norm(epsilon(i,:))^2;%整体误差
end  
end

