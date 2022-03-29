clear
pngname='people2';
RGB=imread([pngname,'.png']); %读取图像
I = rgb2gray(RGB);%转换灰度
BW = imbinarize(I);%二值化
[B,W] = bwboundaries(BW,'holes');%提取边界
imshow(label2rgb(W, @jet, [.5 .5 .5]))%展示连续区域
hold on%展示边界
for k = 1:length(B)
   boundary{k} = B{k};   
   csl = [0;cumsum(vecnorm(diff(boundary{k}),2,2))];
   L(k)=csl(end);
   t{k}=csl/L(k);
   plot(boundary{k}(:,2), boundary{k}(:,1), 'w', 'LineWidth', 2)
end
LL=cumsum(L);
nprate=LL/LL(end);%每段边界长度占总长比例的累加值 用于计算每段取点个数

leastn=ceil(log2(length(B)));%2^leastn>=边界个数
n=7;%预设取点区间为2^n个
N=max([leastn n]);
for i=1:length(nprate)
    newnprate(i)=dec2xiaoshu(xiaoshu2dec(nprate(i),N));
end
while length(newnprate)-length(unique(newnprate))
    N=N+1;
    for i=1:length(nprate)
        newnprate(i)=dec2xiaoshu(xiaoshu2dec(nprate(i),N));
    end
    
end
pointnum=diff([0 newnprate*2^N]);
% while min(pointnum)<=2
%     N=N+1;
%     for i=1:length(nprate)
%         newnprate(i)=dec2xiaoshu(xiaoshu2dec(nprate(i),N));
%     end
%     pointnum=diff([0 newnprate*2^N]);
% end



num=2^N*8;%总点数
figure
hold on
for k = 1:length(B)
    pnum=pointnum(k)*num/2^N;
    %    tt=linspace(0,1,pnum)';
    %     tt=linspace(1,2*pnum-1,pnum)/(pnum*2)';
    %     P{k}=interp1(t{k},boundary{k},tt,'spline');
    %     jianju{k}=vecnorm(diff(P{k}),2,2);
    %     canshucha{k}=diff(tt);
    tt=linspace(0,1,100*length(boundary{k}))';
    PP=interp1(t{k},boundary{k},tt,'makima');
    cslPP=[0;cumsum(vecnorm(diff(PP),2,2))];
    huchang=linspace(1,2*pnum-1,pnum)/(pnum*2)'*cslPP(end);
    for i=1:length(huchang)
        [pMin,minIdx(i)]=min(abs(cslPP-huchang(i)));
    end
    P{k}=PP(minIdx,:);
    
    hcjianju{k}=diff(huchang);
    
    jianju{k}=vecnorm(diff(P{k}),2,2);
    jianju2{k}=diff(cslPP(minIdx));
    clear minIdx
    
    
    scatter(P{k}(:,1),P{k}(:,2));
end
axis equal

pathname='C:\CodeStore\matlab\data\';
Pname=[pngname,num2str(num),'.mat'];

save ([pathname,Pname],'P')



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