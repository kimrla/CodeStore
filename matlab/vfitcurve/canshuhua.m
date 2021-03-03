function t = canshuhua(M,d)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
t(1)=0;%开始参数化 ti i=1,2...,M+1
ds=0;
for j=1:M-1
    ds=ds+norm(d(j+1,:)-d(j,:))^0.5;
end
for i=1:M-1
    t(i+1)=t(i)+norm(d(i+1,:)-d(i,:))^0.5/ds;
end%参数化完成
end

