%% 用V系统特征提取
close all
clear all
example=13;
switch example
    case 1
        % 实验1
        x=0:0.005:1-0.005;
        f=90./(1+exp(-100*(x-0.5)));
        %         M=length(x);%采样点个数
        %         N=32;%基函数个数
        xmax=1;
    case 2
        % 实验2
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-5))+(10*x-5).^5/500;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        xmax=1;
    case 12
        % 实验2 变体3
        x=0:0.02:1-0.02;
        f=100./exp(abs(10*x-5))+(10*x-5).^5/30;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        xmax=1;
    case 11
        % 实验2 变体2
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-5));
        xmax=1;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
    case 13
        % 实验2 变体2
        x=0:0.0002:1-0.0002;
        f=100./exp(abs((10*x-3.5).*(10*x-6.5)));
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        xmax=1;
    case 10
        % 实验2 变体1
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-6))+(10*x-6).^5/500;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        xmax=1;
    case 3
        %
        x=0:0.005:1-0.005;
        f=1./(0.01+(x-0.3).^2).*(x<0.5)+1./(0.015+(x-0.65).^2).*(x>=0.5);
        %         M=length(x);
        %         N=32;%基函数个数
        xmax=1;
    case 4
        x=0:0.02:10-0.02;
        y=@(x) 100./exp(abs(x-5))+(x-5).^5/500;
        wide=1;
        for i=1:length(x)
            f(i)=1/wide*integral(y,x(i)-wide/2,x(i)+wide/2);
        end
        xmax=10;
    case 5
        x=-2:0.04:2-0.04;
        f=sin(x)+2*exp(-30*x.^2);
        x=x+2;
        xmax=4;
    case 6
        x=-2:0.04:2-0.04;
        f=sin(2*x)+2*exp(-16*x.^2)+2;
        x=x+2;
        xmax=4;
    case 7 %间断特征
        x=0:0.005:1-0.005;
        f=4*x.^2.*(3-4*x).*(0<=x&x<0.5)+(4/3*x.*(4*x.^2-10*x+7)-3/2).*(0.5<=x&x<0.75)+(16/3*x.*(x-1).^2).*(0.75<=x&x<=1);
        xmax=1;
    case 8
        x=0:0.005:1-0.005;
        f=2*sin(4*pi*x)-6*(abs(x-0.4)).^0.3-0.5*sign(0.7-x);
        xmax=1;
    case 9
        x=0:0.005:1-0.005;
        f=sin(4*x-2)+2*exp(-30*(4*x-2).^2);
        xmax=1;
    case 14
        load fire500.mat        
end
P=[x;f]';
k=3;
N=floor(log2(length(P)/(k+1)))+1;
t=x';
Lambda=LSMatrix_V(k,N,t)\P;
LambdaA=Lambda;
LambdaA(1:length(Lambda)/2,:)=0;
xijie=LSMatrix_V(k,N,t)*LambdaA;
figure
plot(f)
figure
plot(xijie(:,2))