clear
im=imread('d:\Users\J\图片\实验\shuye1.png');
% imshow(im);
bw = im2bw( im );     

% figure  
% imshow(bw);
% contour = bwperim(bw);
% [m,n]=size(bw);
% bw=ones(m,n)-bw;
% contour = edge(bw ,'Canny');
% % figure
% % imshow(contour)
% contour=rot90(contour,3);
% % contour=ones(m,n)-contour;
B = bwboundaries(bw, 'noholes');
%尺度规范化
[a,~]=size(B);


newB=B{2};
% for i=2:a
%     newB=[newB;B{i}];
% end
maxb=max(max(newB));
% for i=1:a
%     B{i}=B{i}/maxb;
% end
newB=newB/maxb;
% for i=1:a
%     plot(B{i}(:,1),B{i}(:,2));
%     hold on;
% end
% [x,y]=find(contour);
% max=max(max([x,y]));
% x=x/max;
% y=y/max;
% figure
% plot(x,y,'.')
% newB=rot90(newB,3);
% newB=newB';
% figure
% plot(newB(:,1),newB(:,2),'.')
figure
M=200;
for i=1:M
    
    p(i,:)=newB(1+ floor(length(newB)/M*(i-1)),:);    
%     pause(0.1)
%     plot(p(:,1),p(:,2),'.')
end
    x=p(:,2);
    y=1-p(:,1);
    plot(x,y,'.')
    d=[x,y];
    t=canshuhua();


