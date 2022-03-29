function s = SNR(x,y)
%SNR 此处显示有关此函数的摘要
%   此处显示详细说明
% s=10*log10(sum(log10(y.^2./(x-y).^2)));%x原始信号，y重构信号
s=10*log10(sum(y.^2)/sum((x-y).^2));%x原始信号，y重构信号
end

