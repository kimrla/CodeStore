function s = SNR(x,y)
%SNR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% s=10*log10(sum(log10(y.^2./(x-y).^2)));%xԭʼ�źţ�y�ع��ź�
s=10*log10(sum(y.^2)/sum((x-y).^2));%xԭʼ�źţ�y�ع��ź�
end

