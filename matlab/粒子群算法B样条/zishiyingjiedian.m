function ui = zishiyingjiedian(n,p)
%ZISHIYINGJIEDIAN 此处显示有关此函数的摘要
%   此处显示详细说明
ui = zeros(1, n+p+2);
ui(p+2:n+1)=sort(rand(1,n-p));
ui(n+2 : n+p+2) = 1;
end

