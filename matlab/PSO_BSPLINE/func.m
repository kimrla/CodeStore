function result = func(M,n,p,ui,t,l,d)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[N,R,Pi] = kongzhidingdian(M,n,p,t,ui,d);
[epsilon,e] = shujudianwucha(M,N,Pi,d);
[sst,es,ev] = faxiangliangwucha(M,n,p,ui,t,l,Pi);
result=e+ev;
end

