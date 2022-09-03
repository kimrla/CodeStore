% 不提取特征，直接通过加密采样和堆叠基函数的算法
close
clear
plan=2;
switch plan    
    case 1
        load bird200.mat
        N=4;
        rd=3;
    case 2
        load fire500.mat
        N=10;
        rd=2;
        gpoint(end+1,:)=gpoint(1,:);
    case 3
        load yezi600.mat
        N=4;
        rd=1;
        gpoint(end+1,:)=gpoint(1,:);
    case 4
        load shizi1500.mat
        N=4;
        rd=0.1;
    case 5
        load fenghuang2000.mat
        N=10;
%         N=3;
        rd=0.01;
    case 6
        load hudie3000.mat
%         N=10;
        N=3;
        rd=0.01;
    case 7
        load niao21000.mat
%         N=7;
        N=3;
        rd=0.2;
    case 8
        load huacao4-1500.mat
        
%         N=10;
        N=3;
        rd=0.1;
    case 9
        load G-200.mat
%         N=6;
        N=2;
        rd=1;
        gpoint(end+1,:)=gpoint(1,:);
end
P=gpoint;
num=length(P);
k=3;
t=canshuhua(P);

if (k+1)*2^(N-2)+4>num
    tt=unique([linspace(0,1,(k+1)*2^(N-2)+4)';t]);
    P=interp1(t,P,tt,'linear');
    t=tt;
end
VMat=LSMatrix_V(k,N,t);
Lambda=[VMat(:,1:(k+1)) VMat(:,(k+1)*2^(N-2)+1:end)]\P;
curve=[VMat(:,1:(k+1)) VMat(:,(k+1)*2^(N-2)+1:end)]*Lambda;
plot(gpoint(:,1),gpoint(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15,'linewidth',1.5);hold on
plot(curve(:,1),curve(:,2),'Color',[0 102 153]/255,'LineWidth',3)








function [phi,K]=qulv(gpoint,i,num)
r=i+1;
l=i-1;
if i==1
    l=num-1;
elseif i==num
    r=2;
end
v1=gpoint(i,:)-gpoint(l,:);
v2=gpoint(r,:)-gpoint(i,:);
v3=gpoint(r,:)-gpoint(l,:);
phi=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
K=2*sin(phi)/norm(v3);
end

function y=xiaoshu2dec(x,N)
y=zeros(1,N+1);
tempx=x;
for i=1:N+1
    
    y(i)=floor(tempx*2);
    tempx=tempx*2-y(i);
end
if y(end)==1
    y(N)=y(N)+1;
    y(end)=[];
end
for i=fliplr(2:N)
    if y(i)==2
        y(i)=0;
        y(i-1)=y(i-1)+1;
    end
end
% if y(1)==2
%     for i=1:N
%         y(i)=floor(x*2);
%         x=x*2-y(i);
%     end
% end
end

function y=dec2xiaoshu(x)
for i=1:length(x)
    temp(i)=x(i)*2^(-i);
end
y=sum(temp);
end
function [kappa,norm_k] = PJcurvature(gpoint,i)
    x = gpoint(i-1:i+1,1);
    y = gpoint(i-1:i+1,2);
    t_a = norm([x(2)-x(1),y(2)-y(1)]);
    t_b = norm([x(3)-x(2),y(3)-y(2)]);
    
    M =[[1, -t_a, t_a^2];
        [1, 0,    0    ];
        [1,  t_b, t_b^2]];

    a = M\x;
    b = M\y;

    kappa  = 2.*abs(a(3)*b(2)-b(3)*a(2)) / (a(2)^2.+b(2)^2.)^(1.5);
    norm_k =  [b(2),-a(2)]/sqrt(a(2)^2.+b(2)^2.);
end