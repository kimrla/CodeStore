function CC = XGD(x,y)%%��ض�
%XGD �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
avx=mean(x);
avy=mean(y);
CC=sum((x-avx).*(y-avy))/sqrt(sum((x-avx).^2)*sum((y-avy).^2));
end

