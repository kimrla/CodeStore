function CC = XGD(x,y)%%相关度
%XGD 此处显示有关此函数的摘要
%   此处显示详细说明
avx=mean(x);
avy=mean(y);
CC=sum((x-avx).*(y-avy))/sqrt(sum((x-avx).^2)*sum((y-avy).^2));
end

