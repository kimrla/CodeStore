function u= jiedianxiangliang(n,p)
%JIEDIANXIANGLIANG �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
u = zeros(1, n+p+2);%����Ϊ���Ƚڵ�����
piecewise = n - p + 1;       % ���ߵĶ���
if piecewise == 1       % ֻ��һ������ʱ��n = p
    for i = n+2 : n+p+2
        u(i) = 1;
    end
else
    flag = 1;       % ��ֹһ������ʱ
    while flag ~= piecewise
        u(p+1+flag) = u(p + flag) + 1/piecewise;
        flag = flag + 1;
    end
    u(n+2 : n+p+2) = 1;
end

