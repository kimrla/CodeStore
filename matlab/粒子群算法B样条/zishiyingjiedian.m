function ui = zishiyingjiedian(n,p)
%ZISHIYINGJIEDIAN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
ui = zeros(1, n+p+2);
ui(p+2:n+1)=sort(rand(1,n-p));
ui(n+2 : n+p+2) = 1;
end

