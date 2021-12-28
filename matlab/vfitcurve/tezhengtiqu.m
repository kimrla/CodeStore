clear
plan=21;
switch plan
    case 1
        load bird200.mat
        N=4;
        rd=3;
    case 2
        load fire500.mat
        N=4;
        rd=2;
    case 3
        load yezi600.mat
        N=4;
        rd=1;
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
        N=7;
        rd=1;
        
    case 11
        load bird200r.mat
        N=4;
    case 12
        load fire500r.mat
        N=4;
    case 13
        load yezi600r.mat
        N=4;
    case 14
        load shizi1500r.mat
        N=4;
    case 15
        load fenghuang2000.mat
        N=10;
%         N=3;
    case 16
        load hudie3000r.mat
%         N=10;
        N=3;
    case 17
        load niao21000r.mat
%         N=7;
        N=3;
    case 18
        load huacao4-1500r.mat
%         N=10;
        N=3;
    case 19
        load G-200r.mat
%         N=6;
        N=6;
    case 21
        load hudie2fjy420.mat
        N=4;
        gpoint=P;
    case 22
        load star3fjy360.mat
        N=4;
        gpoint=P;
    case 23
        load star3fjy72.mat
        N=4;
        gpoint=P;
end
% gpoint=gpoint+rd*rand(size(gpoint));
% load(['point',num2str(plan),'-200','.mat'])
gpoint(end+1,:)=gpoint(1,:);

t=canshuhua(gpoint);

num=length(gpoint);
for i=1:num
    [phi(i),K(i)]=qulv(gpoint,i,num);
end
q=zeros(num,1);
for i=1:num
    r=i+1;
    l=i-1;
    if i==1
        l=num-1;
    elseif i==num
        r=2;
    end
    if K(i)>K(r) && K(i)>K(l) && phi(i)>pi/6
        %         && phi(i)>pi/6        &&K(i)>3*mean(K)
        q(i)=1;
    end
end

k=3;

tezhengt=t(q==1);

leastduanshu=ceil(log2(length(tezhengt)));%2^leastduanshu>=特征点个数
if leastduanshu<N-1
    leastduanshu=N-1;
end

for i=1:length(tezhengt)
    newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),leastduanshu));
end
while length(newtezhengt)-length(unique(newtezhengt))
    leastduanshu=leastduanshu+1;
    for i=1:length(tezhengt)
        newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),leastduanshu));
    end
end



leastn=leastduanshu+1;
if leastn>N
    N=leastn;
end



% newtezhengt=0:floor((2^(leastduanshu))/(length(tezhengt)-1))/(2^(leastduanshu)):(length(tezhengt)-1)*floor((2^(leastduanshu))/(length(tezhengt)-1))/(2^(leastduanshu));

% newtezhengt=0:floor((2^(N-1))/(length(tezhengt)-1))/(2^(N-1)):(length(tezhengt)-1)*floor((2^(N-1))/(length(tezhengt)-1))/(2^(N-1));


newt(q==1)=newtezhengt;
tlist=find(q);

%
% newt(fenduandian)=newtezhengt;
% tlist=fenduandian;

for i=2:length(tlist)
    j=tlist(i-1)+1:tlist(i)-1;
    newt(j)=(t(j)-t(tlist(i-1)))*(newtezhengt(i)-newtezhengt(i-1))/(tezhengt(i)-tezhengt(i-1))+newt(tlist(i-1));
end
% newt(tlist(end)+1:num)=(t(tlist(end)+1:num)-t(tlist(end)))*(1-newtezhengt(end))/(1-tezhengt(end))+newt(tlist(end));

newt=newt';

CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));

CList(:,2)=CList(:,1);
CList(:,3)=2;
CList(ismember(CList(:,1),newtezhengt),3)=0;
find(CList(:,3)==0)
CList(end+1,:)=[0,1,0];
Lambda = LSCurFit_V(gpoint,k,N,newt,CList);
Lambda(all(abs(Lambda)<=10^(-3),2),:)=0;
nc=LSMatrix_V(k,N,newt)*Lambda;
% wucha=vecnorm((nc-gpoint),2,2);
% pjwucha=mean(wucha);
[wuchaV,pjwuchaV]=distanceerror(gpoint,nc);


pathname='C:\CodeStore\matlab\vfitcurve\data\';

tlistname=['tlist',num2str(plan),'.mat'];
save ([pathname,tlistname],'tlist','plan')
tzwuchaV=wuchaV(tlist);
tzwcname=['tzwuchaV',num2str(plan),'.mat'];
save ([pathname,tzwcname],'tzwuchaV','pjwuchaV','wuchaV') 



figure,
if plan<20
load(['point',num2str(plan),'-200','.mat'])
end
plot(gpoint(:,1),gpoint(:,2),'.','Color','r','MarkerSize',10);hold on
% plot(gpoint(tlist(21:end),1),gpoint(tlist(21:end),2),'.','Color',[224 222 58]/255,'MarkerSize',15);hold on

% gpoint(end,:)=NaN;
% patch(gpoint(:,1),gpoint(:,2),K','EdgeColor','interp','Marker','.','MarkerSize',15,'MarkerFaceColor','flat');hold on
% colormap(magma)
% 
% fenduan=1/(2^(N-2)):1/(2^(N-2)):(2^(N-2)-1)/(2^(N-2));
% A=LSMatrix_V(k,N,fenduan');
% jiedian=A*Lambda;
% T=LSMatrix_V(k,N,newtezhengt')*Lambda;

% scatter(jiedian(:,1),jiedian(:,2),'b')
% hold on
% plot(T(:,1),T(:,2),'gs','color',[0 102 153]/255,'MarkerSize',10)
VCompose(Lambda,k,N)
% legend('原始数据点','特征点','拟合曲线')
legend({'原始数据','拟合曲线'},'location','northwest','fontsize', 15, 'fontname', '微软雅黑')
% set(gca, 'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
axis equal
axis off
% plot(gpoint(u,1),gpoint(u,2),'.','Color',[255 0 102]/255,'MarkerSize',15);
Lambda(all(Lambda~=0,2),3)=find(all(Lambda~=0,2));
Lambda(all(Lambda==0,2),:)=[];


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

% function [a,b]=jianruijida(gpoint,k,t)
%     N=floor(log2(length(gpoint)));
%     A=LSMatrix(k,N,t);
%     Lambda=A\gpoint;
%     DC3=Lambda((k+3)*2^(N-2)+1:(k+4)*2^(N-2),:);
%     [~,locationx] = findpeaks(abs(DC3(:,1)));
%     [~,locationy] = findpeaks(abs(DC3(:,2)));
%     xa=(locationx-1)/2^(N-2);
%     xb=locationx/2^(N-2);
%     ya=(locationy-1)/2^(N-2);
%     yb=locationy/2^(N-2);
%     a=[xa,ya];
%     b=[xb,yb];
% end

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

function Kse=qulvjinsi(K,s,e,u)
i=s:e-1;
Kse=sum((K(i)+K(i+1)).*(u(i+1)-u(i))'/2);
end

function Lse=huchangjinsi(s,e,gpoint)
Lse=sum(vecnorm(diff(gpoint(s:e)),2,2));
end

function SIse=SIse(K,s,e,u,gpoint)
r=0.5;
SIse=r*qulvjinsi(K,s,e,u)/qulvjinsi(K,1,length(gpoint),u)+(1-r)*huchangjinsi(s,e,gpoint)/huchangjinsi(1,length(gpoint),gpoint);
end

