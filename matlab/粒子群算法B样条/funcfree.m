function result = funcfree(M,n,p,ui,t,d)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[N,R,Pi] = kongzhidingdian(M,n,p,t,ui,d);
[epsilon,e] = shujudianwucha(M,N,Pi,d);
result=(M+1)*log(e)+(log(M+1))*(2*n-p+1);
end

