clear
load yezi600.mat
t=huancanshuhua(gpoint);
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
    if K(i)>K(r) && K(i)>K(l) && phi(i)>pi/4
        q(i)=1;
    end
end
tezhengt=t(q==1);
leastduanshu=ceil(log2(length(tezhengt)));

for i=1:length(tezhengt)
    newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),leastduanshu));
end
while length(newtezhengt)-length(unique(newtezhengt))
    leastduanshu=leastduanshu+1;
    for i=1:length(tezhengt)
    newtezhengt(i)=dec2xiaoshu(xiaoshu2dec(tezhengt(i),leastduanshu));
    end
end
newt(q==1)=newtezhengt;
tlist=find(q);
for i=2:length(tlist)
    j=tlist(i-1)+1:tlist(i)-1;
    newt(j)=(t(j)-t(tlist(i-1)))*(newtezhengt(i)-newtezhengt(i-1))/(tezhengt(i)-tezhengt(i-1))+newt(tlist(i-1));   
end
    newt(tlist(end)+1:num)=(t(tlist(end)+1:num)-t(tlist(end)))*(1-newtezhengt(end))/(1-tezhengt(end))+newt(tlist(end));
    leastn=leastduanshu+1;

newt=newt';
k=3;
N=5;
if leastn>N
    N=leastn;
end

% A = LSMatrix_V(k,N,newt);
[DR,DL] = VContinuityInfo1(N);
CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=2;
CList(ismember(CList(:,1),newtezhengt),3)=0;
Lambda = LSCurFit_V(gpoint,k,N,newt,CList);
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
%         plot(nx,ny,'.')


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
VCompose(Lambda,k,N)
axis equal
axis off













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
