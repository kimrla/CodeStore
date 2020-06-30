function u = jiedianxiangliang3(n,p,t,M)
%JIEDIANXIANGLIANG3 此处显示有关此函数的摘要
%   此处显示详细说明
u = zeros(1, n+p+2);
for j=1:n-p
    a=(M+1)/(n-p+1);
    i=floor(j*a);
    u(j+p+1)=(1-(j*a-i))*t(i)+(j*a-i)*t(i+1);
end
u(n+2:n+2+p)=1;
end

