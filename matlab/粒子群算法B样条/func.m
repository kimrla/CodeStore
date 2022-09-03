function result = func(M,n,p,ui,t,l,d)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
[N,R,Pi] = kongzhidingdian(M,n,p,t,ui,d);
[epsilon,e] = shujudianwucha(M,N,Pi,d);
[sst,es,ev] = faxiangliangwucha(M,n,p,ui,t,l,Pi);
result=e+ev;
end

