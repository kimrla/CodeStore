clear
RGB=imread('people2.png'); %读取图像
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
beilv=L/min(L);
decbeilv=2.^round(log2(beilv));


num=1000;%总点数
figure
hold on
for k = 1:length(B)
   pointnum=round(num*nprate(k));
   tt=linspace(0,1,pointnum)';  
   P{k}=interp1(t{k},boundary{k},tt,'linear');
   scatter(P{k}(:,1),P{k}(:,2));
end
axis equal

