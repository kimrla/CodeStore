function u= jiedianxiangliang2(n,p,t,M)
%JIEDIANXIANGLIANG2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
u = zeros(1, n+p+2);
for j=1:n-p
        a=(M+1)/(n-p+1);
        i=floor(j*a)+1;
%     for i=j:j+p-1
        u(j+p+1)=u(j+p+1)+1/p*sum(t(i:i+p-1));
%     end
end
u(n+2:n+2+p)=1;
end

