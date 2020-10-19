function [BIC] = BIC(Num,M,n,p,t,ui,d)
%BIC 此处显示有关此函数的摘要
%   此处显示详细说明
[N,~,P] = kongzhidingdian(M,n,p,t,ui,d);
[~,R]=shujudianwucha(M,N,P,d);
BIC=Num*log(R)+log(Num)*(p+n+1); 
end

