
%% ��1:��һ�ι⻬����(����)��������1��V����С���˱ƽ�

clear all
close all

%% Test1
x = linspace(0,1,100)';
Y = cos(4*pi*x);
figure,
plot(x,Y);

k = 1;
N = 5;
A = LSMatrix_V(k,N,x);

X = A\Y;

V = LSMatrix_V(k,N,linspace(0,1,501)');
F = V * X;
figure,
plot(x,Y,'.r');hold on
plot(linspace(0,1,501)',F);