function [BIC] = BIC(Num,M,n,p,t,ui,d)
%BIC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[N,~,P] = kongzhidingdian(M,n,p,t,ui,d);
[~,R]=shujudianwucha(M,N,P,d);
BIC=Num*log(R)+log(Num)*(p+n+1); 
end

