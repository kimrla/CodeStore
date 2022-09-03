function result = funcfree(M,n,p,ui,t,d)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
[N,R,Pi] = kongzhidingdian(M,n,p,t,ui,d);
[epsilon,e] = shujudianwucha(M,N,Pi,d);
result=(M+1)*log(e)+(log(M+1))*(2*n-p+1);
end

