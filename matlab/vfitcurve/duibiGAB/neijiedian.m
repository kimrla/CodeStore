function u= neijiedian(n,p,a,b)
%JIEDIANXIANGLIANG 此处显示有关此函数的摘要
%   此处显示详细说明
u = zeros(1, n-p);%以下为均匀节点向量
piecewise = n - p + 1;       % 曲线的段数
% if piecewise == 1       % 只有一段曲线时，n = p
%     for i = n+2 : n+p+2
%         u(i) = 1;
%     end
% else
    u(1)=a+b/piecewise;
    flag = 1;       % 不止一段曲线时
    while flag ~= piecewise-1
        u(1+flag) = u(flag) + b/piecewise;
        flag = flag + 1;
    end    
end

