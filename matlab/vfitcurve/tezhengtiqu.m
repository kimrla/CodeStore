clear
load fire200.mat
t=canshuhua(num,gpoint);
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
    if K(i)>K(r) && K(i)>K(l) && phi(i)>pi/2
        q(i)=1;
    end
end
tezhengt=t(q==1);






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
