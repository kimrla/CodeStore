function [sst,es,ev] = faxiangliangwucha(M,n,p,ui,t,l,P)
%FAXIANGLIANGWUCHA 此处显示有关此函数的摘要
%   此处显示详细说明
sst=zeros(M+1,2);
for i=1:M+1
    for j=1:n+1
        Length1 = ui(j+p) - ui(j);
        Length2 = ui(j+p+1) - ui(j+1);      % 支撑区间的长度
        if Length1 == 0.0       % 规定0/0 = 0
            Length1 = 1.0;
        end
        if Length2 == 0.0
            Length2 = 1.0;
        end
        sst(i,:)=sst(i,:)+p*(Njp(j,p-1,t(i),ui)/Length1-Njp(j+1,p-1,t(i),ui)/Length2)*P(j,:);
    end
end
ev=0;
for i=1:M+1
    es(i)=(dot(sst(i,:),l(i,:)))^2;
    ev=ev+es(i);
end
end

