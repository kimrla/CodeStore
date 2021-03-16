clear
example=1;
switch example
    case 1
        % 实验1
        x=0:0.005:1-0.005;
        f=90./(1+exp(-100*(x-0.5)));        
        M=length(x);%采样点个数
        N=32;%基函数个数
        xmax=1;
    case 2
        % 实验2
        x=0:0.05:10-0.05;
        f=100./exp(abs(x-5))+(x-5).^2/500;        
        M=length(x);%采样点个数  
        N=128;%基函数个数
        xmax=10;
    case 3
        %
        x=0:0.005:1-0.005;
        f=1./(0.01+(x-0.3).^2).*(x<0.5)+1./(0.015+(x-0.65).^2).*(x>=0.5);
        M=length(x);
        N=32;%基函数个数
        xmax=1;
end
f_=f+normrnd(0,1,1,M);

k=1;%次数
% t=canshuhua(M,d);
t=x/xmax;
d=[t;f_]';
[V,lamda]=vxishu(k,N,M,t,d);
[V1,lamda1]=vdxishu(k,N,M,t,d);
plot(t,f,'*')
hold on
drawV(N, k, lamda);
% figure
% for i=1:N
%     subplot(N,1,i),plot(x,V(:,i))
% end


