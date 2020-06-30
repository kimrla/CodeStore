function u= jiedianxiangliang2(n,p,t)
%JIEDIANXIANGLIANG2 此处显示有关此函数的摘要
%   此处显示详细说明
u = zeros(1, n+p+2);
for j=1:n-p
    for i=j:j+p-1
        u(j+p+1)=u(j+p+1)+1/p*t(i);
    end
end
u(n+2:n+2+p)=1;
end

