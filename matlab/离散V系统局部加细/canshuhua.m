function t = canshuhua(d)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
% t(1)=0;%开始参数化 ti i=1,2...,M+1
% ds=0;
% M=length(d);

% for j=1:M-1
%     ds=ds+norm(d(j+1,:)-d(j,:))^0.5;
% end
% % ds=ds+norm(d(1,:)-d(end,:))^0.5;
% % j=1:M-1;
% % ds=sum(norm(d(j+1,:)-d(j,:))^0.5);
% for i=1:M-1
%     t(i+1)=t(i)+norm(d(i+1,:)-d(i,:))^0.5/ds;
% end%参数化完成
csl = [0;cumsum(vecnorm(diff(d),2,2).^0.5)];       % 累加长度
LChord = csl(end);
t = csl / LChord;
end

