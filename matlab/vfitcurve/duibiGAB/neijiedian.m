function u= neijiedian(n,p,a,b)
%JIEDIANXIANGLIANG �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
u = zeros(1, n-p);%����Ϊ���Ƚڵ�����
piecewise = n - p + 1;       % ���ߵĶ���
% if piecewise == 1       % ֻ��һ������ʱ��n = p
%     for i = n+2 : n+p+2
%         u(i) = 1;
%     end
% else
    u(1)=a+b/piecewise;
    flag = 1;       % ��ֹһ������ʱ
    while flag ~= piecewise-1
        u(1+flag) = u(flag) + b/piecewise;
        flag = flag + 1;
    end    
end
