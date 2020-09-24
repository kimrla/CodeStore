clear all;
close all; 
clc;
x=0:1:10; 
y=[-0.447,1.978,3.28,6.16,7.08,7.34,7.66,9.56,9.48,9.3,11.2]; 
plot(x,y,'k.','markersize',25) 
hold on                                       %图形保持功能  
axis([0 13 -2 16])                           %坐标轴范围 
order=3;       %阶数（一般不要超过3，可以更改确定） 
[p, s ]=polyfit(x,y,order);
t=0:1:12; 
Q=polyval(p,t); 
plot(t,Q,'r');                                    %画出图形
r=corrcoef(x,y);                             %相关系数

