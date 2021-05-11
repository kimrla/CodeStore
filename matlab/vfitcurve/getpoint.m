% 绘制三种类型的B样条曲线，需要前面所给的所有.m文件
clear all;
%控制顶点
P{1}=[91.5,414.2;122,244;232,86;358.5,9];
P{2}=[358.5,9;260,172;405,277;357.4,479.5];
P{3}=[357.4,479.5;398,412;428,333;421.1,201.3];
P{4}=[421.1,201.3;482,258;509,324;529.3,406.4];
P{5}=[529.3,406.4;543,386;554,366;567.1,358.3];
P{6}=[567.1,358.3;525,502;728,515;459.9,751.1];
P{7}=[459.9,751.1;544,571.5;475,485;443.9,456.4];
P{8}=[443.9,456.4;425,538.5;393,614.4;341.6,666.6];
P{9}=[341.6,666.6;318,502;269,555;256.2,292.2];
P{10}=[256.2,292.2;234,349;199.5,444;185,493.4];
P{11}=[185,493.4;176,481;163,462;155.1,460.7];
P{12}=[155.1,460.7;195.5,561.5;58.5,618;159,746.2];
P{13}=[159,746.2;-121,530;96,486;46,358.2];
P{14}=[46,358.2;69.5,380;73,384;91.5,414.2];


p = 3;
num=200;

for i=1:length(P)
cp{i,1}=redraw(p,P{i});%根据控制顶点重新绘制曲线
end


curve=cell2mat(cp);
curve(:,2)=max(curve(:,2))-curve(:,2);
for i=1:length(P)
    cp{i}(:,2)=max(curve(:,2))-cp{i}(:,2);
end
% jianju(1)=norm(curve(end,:)-curve(1,:));
% l(1)=jianju(1);
% for i=2:length(curve)
%     jianju(i)=norm(curve(i,:)-curve(i-1,:));
%     l(i)=l(i-1)+jianju(i);
% end
jianju(length(curve))=norm(curve(end,:)-curve(1,:));

for i=2:length(curve)
    jianju(i-1)=norm(curve(i,:)-curve(i-1,:));    
end
l(1)=jianju(1);
for i=2:length(curve)
     l(i)=l(i-1)+jianju(i);  
end

for i=1:length(cp)
%     j=1+(i-1)*200:i*200;
    lcp(1,i)=jianju(1+(i-1)*200);
    for j=2:200
        lcp(j,i)=lcp(j-1,i)+jianju(j+(i-1)*200);
    end
%     j=1+(i-1)*200:i*200;
%     lcp(i)=sum(jianju(j));
    ratecp(i)=lcp(end,i)/l(end);
    gpcp{i,1}=cp{i}(1,:);
end
ncptemp=ratecp*num;
ncp(1)=round(ncptemp(1));
for i=2:length(ncptemp)-1
    ncp(i)=round(ncp(i-1)+ncptemp(i))-ncp(i-1);
end
ncp(end+1)=num-sum(ncp);

cpinterval=lcp(end,:)./ncp;
for i=1:size(lcp,2)
    counter=1;
    for j=1:size(lcp,1)-1
        if floor(lcp(j,i)/(counter*cpinterval(i)))
            counter=counter+1;
%             if j==size(lcp,1)
%             else
            gpcp{i}(counter,:)=cp{i}(j+1,:);    
%             end
        end
    end
end

gpoint=cell2mat(gpcp);
% 
% for i=1:length(gpoint)
%     pause(0.1)
%     plot(gpoint(i,1),gpoint(i,2),'.')
%     hold on
% end









% gpinterval=l(end)/num;
% for i=1:size(lcp,2)
%     counter=1;
%     for j=1:size(lcp,1)
%         if floor(lcp(j,i)/(counter*gpinterval))
%             counter=counter+1;
%             gpcp{i}(counter,:)=cp{i}(j,:);
% %             if lcp(end,i)-lcp(j,i)<gpinterval
% %                 gpcp{i}(counter,:)=[];
% %             end
%         end
%     end
% end
% if norm(gpcp{1}(1,:)-gpcp{end}(end,:))<gpinterval
%     gpcp{end}(end,:)=[];
% end
% for i=2:length(gpcp)
%     if norm(gpcp{i}(1,:)-gpcp{i-1}(end,:))<gpinterval
%         gpcp{i-1}(end,:)=[];
%     end
% end
    


















% for i=1:length(curve)
% pause(0.02)
% plot(curve(i,1),curve(i,2),'.');
% hold on
% end

function p_u=redraw(p,P)

n = length(P)-1;  

ui=jiedianxiangliang(n,p,0,1); % 准均匀B样条的节点矢量

Njp_u = zeros(1, n+1);
j=0;
for u = 0 : 0.005 : 1-0.005
    j=j+1;
%     p_u=zeros(length(u),2);
    for i = 1 : n+1
        Njp_u(j, i) = Njp(i, p , u, ui);
    end
    p_u(j,:) = Njp_u(j,:)*P;       
    
end
end
