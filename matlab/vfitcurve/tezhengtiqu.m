clear
plan=5;
switch plan
    case 1
        load bird200.mat
        N=4;
    case 2
        load fire500.mat
        N=4;
    case 3
        load yezi600.mat
        N=4;
    case 4
        load shizi1500.mat
        N=4;
    case 5
        load fenghuang2000.mat
        N=10;
end
% load zsy.mat
% gpoint=new_C{1};
t=huancanshuhua(gpoint);
% gpoint=[gpoint;gpoint(1,:)];
% t=canshuhua(gpoint);
% gpoint=gpoint/max(max(gpoint));
num=length(gpoint);
for i=1:num
    [phi(i),K(i)]=qulv(gpoint,i,num);
end
q=zeros(num,1);
for i=1:num
    r=i+1;
    l=i-1;
    if i==1
        l=num;
    elseif i==num
        r=1;
    end
    if K(i)>K(r) && K(i)>K(l) && phi(i)>pi/6
        %         && phi(i)>pi/6        &&K(i)>3*mean(K)
        q(i)=1;
    end
end

tezhengt=t(q==1);
% fenduandian(end+1)=50;
% fenduandian(end+1)=2001;
% fenduandian(end+1)=1950;

% fenduandian=sort(fenduandian);
% tezhengt=t(fenduandian);



leastduanshu=ceil(log2(length(tezhengt)));


%
for i=1:length(tezhengt)
    newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),leastduanshu));
end
while length(newtezhengt)-length(unique(newtezhengt))
    leastduanshu=leastduanshu+1;
    for i=1:length(tezhengt)
        newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),leastduanshu));
    end
end

k=3;

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
newt(tlist(end)+1:num)=(t(tlist(end)+1:num)-t(tlist(end)))*(1-newtezhengt(end))/(1-tezhengt(end))+newt(tlist(end));

newt=newt';


% A = LSMatrix_V(k,N,newt);
% [DR,DL] = VContinuityInfo1(N);
CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=2;
CList(ismember(CList(:,1),newtezhengt),3)=0;
CList(end+1,:)=[0,1,0];
Lambda = LSCurFit_V(gpoint,k,N,newt,CList);
Lambda(all(abs(Lambda)<=10^(-3),2),:)=0;
[fitc,cost]=wucha(k,N,newt,Lambda,gpoint);
meanwc=cost(end)/length(gpoint);

tlistname=['tlist',num2str(plan),'.mat'];
save (tlistname,'tlist','plan')

tezhengc=LSMatrix_V(k,N,newtezhengt')*Lambda;

tzwuchaV=vecnorm((tezhengc-gpoint(tlist,:)),2,2);
tzwcname=['tzwuchaV',num2str(plan),'.mat'];
save (tzwcname,'tzwuchaV') 
% for i=1:length(tlist)
tlist(end+1)=length(gpoint);
wcrate=diff(cost(tlist))./diff(tlist);
% end

[~,jc]=max(wcrate);
u=tlist(jc)+1:tlist(jc+1)-1;
% % 
% N=N+1;
% CListp(:,1)=newt(tlist(jc)):1/(2^(N-1)):newt(tlist(jc+1));
% CListp(:,2)=CListp(:,1);
% CListp(:,3)=2;
% FP=fitc(u,:)-gpoint(u,:);
% Lplus=fitjiaceng(FP,k,N,newt,u,CListp);
% 
% Lambda=[Lambda;Lplus];
% [~,costs]=wucha(k,N,newt,Lambda,gpoint);
% wcrate2=diff(costs(tlist))./diff(tlist);

% for i=tlist(jc)+1:tlist(jc+1)-1
% SIsi(i)=abs(SIse(K,tlist(jc),i,newt,gpoint)-0.5*SIse(K,tlist(jc),tlist(jc+1),newt,gpoint));
% end
% [~,Pnew]=min(SIsi(tlist(jc)+1:tlist(jc+1)-1));
% Pnew=Pnew+tlist(jc);
% newt(Pnew)







% p = sum(CList(:,3)+1);  % 约束方程数目
%         SegNum = 2^(N-1); % 分段数
%         VNum = 2^N;       % 基函数数目
%         knots = linspace(0,1,SegNum+1);% 节点向量
%         C = zeros(p,VNum);    % 约束矩阵
%         csidx = 0;
%         for h = 1 : length(CList)           % 给每个约束建立相应的方程组
%             x1 = find(knots==CList(h,1));   % 左节点
%             x2 = find(knots==CList(h,2));   % 右节点
%             for c = 0 : CList(h,3)          % 从C^0约束至C^CList(h,3)约束，每个约束建立一个方程
%                 csidx = csidx + 1;
%                 C(csidx,:) = [DR(:,x1,c+1) - DL(:,x2,c+1)]';
%             end
%         end
%
%         M = [2*A'*A,C';
%             C,zeros(p)];
%         d = zeros(p,1);
%         bx = [2*A'*gpoint(:,1);d];
%         by = [2*A'*gpoint(:,2);d];
%
%         X = M\bx;
%         Y = M\by;
%
%         X = X(1:VNum);
%         Y = Y(1:VNum);
%         figure
%         A = LSMatrix_V(k,N,newt);
%         nx=A*Lambda(:,1);
%         nx(end+1)=nx(1);
%         ny=A*Lambda(:,2);
%         ny(end+1)=ny(1);
%         plot(gpoint(:,1),gpoint(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
%         plot(nx,ny)

% figure
% caiyang=0:1/2000:1;
% plot(gpoint(:,1),gpoint(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
% curve=LSMatrix_V(k,N,caiyang')*Lambda;
% plot(curve(:,1),curve(:,2),'b')
% fenduan=1/(2^(N-2)):1/(2^(N-2)):(2^(N-2)-1)/(2^(N-2));
% A=LSMatrix_V(k,N,fenduan');
% jiedian=A*Lambda;
% T=LSMatrix_V(k,N,newtezhengt')*Lambda;
% scatter(jiedian(:,1),jiedian(:,2),'b')
% % hold on
% scatter(T(:,1),T(:,2),'g')




% figure
% bar(abs(Y))
% % 找到最大的lamda
%
% matrixLamda=ones(length(Y)/2,N);
% matrixLamda(:,1)=kron(Y(1:2),ones(32,1));
% for i=2:N
% matrixLamda(:,i)=kron(Y(1+(k+1)*2^(i-2):(k+1)*2^(i-1)),ones(length(Y)/2^i,1));
% end
% matrixLamda=abs(matrixLamda);
% imagesc(matrixLamda(:,2:N))
% colorbar
%
% InfoV1Bas = BaseGene_V1(N);   % 1次V系统基函数信息(非离散采样)
% [VRInfo,NumSeg] = VReconstruction_Polyline([X Y],InfoV1Bas);
% figure,
% plot(gpoint(:,1),gpoint(:,2),'.','Color',[255 102 102]/255,'MarkerSize',10);hold on
% for j = 1 : NumSeg
%     plot(VRInfo(j,5:6),VRInfo(j,7:8),'Color',[0 102 153]/255,'LineWidth',1.5);hold on
% end
% nx=A*X;
% nx(end+1)=nx(1);
% ny=A*Y;
% ny(end+1)=ny(1);
% plot(nx,ny,'Color',[0 102 153]/255,'LineWidth',1.5)


figure,
plot(gpoint(:,1),gpoint(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
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
legend('原始数据点','拟合曲线')
axis equal
axis off
% plot(gpoint(u,1),gpoint(u,2),'.','Color',[255 0 102]/255,'MarkerSize',15);
Lambda(all(Lambda~=0,2),3)=find(all(Lambda~=0,2));
Lambda(all(Lambda==0,2),:)=[];

% nodes = [0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
% treeplot(nodes,'[0 0 1]')










function [phi,K]=qulv(gpoint,i,num)
r=i+1;
l=i-1;
if i==1
    l=num;
elseif i==num
    r=1;
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
if y(1)==2
    for i=1:N
        y(i)=floor(x*2);
        x=x*2-y(i);
    end
end
end

function y=dec2xiaoshu(x)
for i=1:length(x)
    temp(i)=x(i)*2^(-i);
end
y=sum(temp);
end

function [c,y]=wucha(k,N,t,lambda,P)
c=LSMatrix_V(k,N,t)*lambda;
y=cumsum(vecnorm((c-P),2,2));
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

% im=imread('d:\Users\J\图片\实验\xiaoniao1.png');
% imshow(im);
% bw = im2bw( im );
%
% % figure
% % imshow(bw);
% % contour = bwperim(bw);
% % [m,n]=size(bw);
% % bw=ones(m,n)-bw;
% % contour = edge(bw ,'Canny');
% % % figure
% % % imshow(contour)
% % contour=rot90(contour,3);
% % % contour=ones(m,n)-contour;
% B = bwboundaries(bw, 'noholes');
% %尺度规范化
% [a,~]=size(B);
%
%
% newB=B{2};
% % for i=2:a
% %     newB=[newB;B{i}];
% % end
% maxb=max(max(newB));
% % for i=1:a
% %     B{i}=B{i}/maxb;
% % end
% newB=newB/maxb;
% % for i=1:a
% %     plot(B{i}(:,1),B{i}(:,2));
% %     hold on;
% % end
% % [x,y]=find(contour);
% % max=max(max([x,y]));
% % x=x/max;
% % y=y/max;
% % figure
% % plot(x,y,'.')
% % newB=rot90(newB,3);
% % newB=newB';
% figure
% plot(newB(:,1),newB(:,2),'.')
% figure
% M=200;
% for i=1:M
%
%     p(i,:)=newB(1+ floor(length(newB)/M*(i-1)),:);
% %     pause(0.1)
% %     plot(p(:,1),p(:,2),'.')
% end
%     x=p(:,2);
%     y=1-p(:,1);
%     plot(x,y,'.')
% %     d=[x,y];
% %     d=d';
% %     t=canshuhua(M,d);
% %
% %     k=1;
% %     N=7;
% %
% %
% %
% %
% %
% %
