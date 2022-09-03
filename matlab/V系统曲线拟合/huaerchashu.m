%% 画LAMBDA系数的二叉树和热力图
clear
close all
example=10;
switch example
    case 1
        % 实验1
        x=0:0.005:1-0.005;
        f=90./(1+exp(-100*(x-0.5)));
        %         M=length(x);%采样点个数
        %         N=32;%基函数个数
        
    case 2
        % 实验2
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-5))+(10*x-5).^5/500;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
       
    case 12
        % 实验2 变体3
        x=0:0.02:1-0.02;
        f=100./exp(abs(10*x-5))+(10*x-5).^5/30;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        
    case 11
        % 实验2 变体2
        x=0:0.002:1-0.002;
        f=100./exp(abs(10*x-5));
        
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
    case 13
        % 实验2 变体2
        x=0:0.0002:1-0.0002;
        f=100./exp(abs((10*x-2.5).*(10*x-7.5)));
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
        
    case 10
        % 实验2 变体1
        x=0:0.002:1;
        f=100./exp(abs(10*x-6.5))+(10*x-6.5).^5/500;
        %         M=length(x);%采样点个数
        %         N=128;%基函数个数
      
    case 3
        %
        x=0:0.005:1-0.005;
        f=1./(0.01+(x-0.3).^2).*(x<0.5)+1./(0.015+(x-0.65).^2).*(x>=0.5);
        %         M=length(x);
        %         N=32;%基函数个数
        
    case 4
        x=0:0.02:10-0.02;
        y=@(x) 100./exp(abs(x-5))+(x-5).^5/500;
        wide=1;
        for i=1:length(x)
            f(i)=1/wide*integral(y,x(i)-wide/2,x(i)+wide/2);
        end
       
    case 5
        x=-2:0.04:2-0.04;
        f=sin(x)+2*exp(-30*x.^2);
        x=x+2;
      
    case 6
        x=-2:0.04:2-0.04;
        f=sin(2*x)+2*exp(-16*x.^2)+2;
        x=x+2;
        
    case 7 %间断特征
        x=0:0.005:1-0.005;
        f=4*x.^2.*(3-4*x).*(0<=x&x<0.5)+(4/3*x.*(4*x.^2-10*x+7)-3/2).*(0.5<=x&x<0.75)+(16/3*x.*(x-1).^2).*(0.75<=x&x<=1);
       
    case 8
        x=0:0.005:1-0.005;
        f=2*sin(4*pi*x)-6*(abs(x-0.4)).^0.3-0.5*sign(0.7-x);
     
    case 9
        x=0:0.005:1-0.005;
        f=sin(4*x-2)+2*exp(-30*(4*x-2).^2);
      
end
plot(x,f,'.','Color',[255 102 102]/255,'MarkerSize',15)
%% 普通V系统拟合（非正交变换，无特征提取，无约束条件，不参数化）
k=3;
N=5;

lambda=LSMatrix_V(k,N,x')\f';
figure
plot(x,f,'.','Color',[255 102 102]/255,'MarkerSize',15),hold on;
VCompose(lambda,k,N)

%% 画系数lambda的树状图
s = 1:k+1;
t = k+2:(k+1)*2;
if N>2
    s=[s,kron(k+2:length(lambda)/2,ones(1,2))];
    t=[t,(k+1)*2+1:length(lambda)];
end

% normlambda(1:k+1)=abs(lambda(1:k+1))/sum(abs(lambda(1:k+1)));
% normlambda(k+2:(k+1)*2)=abs(lambda(k+2:(k+1)*2))/sum(abs(lambda(k+2:(k+1)*2)));
normlambda(1:(k+1)*2)=1;
intervalsize(1:(k+1)*2)=1;
for n=3:N %将每一类基函数的系数分别取绝对值后归一化
    for j=1:k+1
    normlambda((k+j)*2^(n-2)+1:(k+1+j)*2^(n-2))= ...,
        abs(lambda((k+j)*2^(n-2)+1:(k+1+j)*2^(n-2))) ...,
        /sum(abs(lambda((k+j)*2^(n-2)+1:(k+1+j)*2^(n-2))));
    end
    intervalsize((k+1)*2^(n-2)+1:(k+1)*2^(n-1))=1/2^(n-2);
end

G = digraph(s,t);
figure
tree=plot(G,'MarkerSize',4,'Layout','Layered');
set(tree,'NodeCData',normlambda)
set(tree,'MarkerSize',intervalsize*110)
tree.Marker='s';
colormap()
colorbar
%% 热力图
heatlambda=zeros(length(lambda)/2,N);
normheatlambda=heatlambda;
heatlambda(:,1:2)=kron(reshape(lambda(1:(k+1)*2),[],2),ones(2^(N-2),1));
for n=3:N
    heatlambda(:,n)=kron(lambda((k+1)*2^(n-2)+1:(k+1)*2^(n-1)),ones(2^(N-n),1));
end
heatlambda=abs(heatlambda);
normheatlambda(:,1:2)=1;
for n=3:N  %将每一类基函数的系数分别取绝对值后归一化
    for j=1:k+1
        normheatlambda(1+(j-1)*2^(N-2):j*2^(N-2),n)= ...,
            heatlambda(1+(j-1)*2^(N-2):j*2^(N-2),n)/ ...,
            sum(abs(lambda((k+j)*2^(n-2)+1:(k+1+j)*2^(n-2))));
    end
end
figure
for j=1:k+1
ylable=[0 1];
xlable=1:N;
subplot(k+1,1,j),imagesc(xlable,ylable,normheatlambda(1+(j-1)*2^(N-2):j*2^(N-2),:))
set(gca,'xtick',1:N)
end
