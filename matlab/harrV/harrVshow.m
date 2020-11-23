% x1=linspace(0,0.5);
% x2=linspace(0.5,1);
% y1=x1;
% y2=-x2;
% plot(x1,y1,x2,y2)
clear all
hold on
% plot(x1,y1)
% plot(x2,y2)
x0=linspace(0,1);
y0=8^0.5*3^0.5*(1-2*x0);
% y0=8^0.5*5^0.5*(6*x0.^2-6*x0+1);
% y0=8^0.5*7^0.5*(-20*x0.^3+30*x0.^2-12*x0+1);
% plot(x0,y0);
for i=1:8
    x{i}=linspace((i-1)/8,i/8);
    y{i}=y0;   
    plot(x{i},y{i},'b')
end
for i=1:4
    P{i}=hadamard(2^(i-1));
end
% for i=1:4
%     T{i}=zeros(8,8);
% end
for j=1:4
    T{j}=P{j};
for i=1:8/2^(j-1)-1
    T{j}=blkdiag(T{j},P{j});
end
end
        % T{2}=blkdiag(P{2},P{2},P{2},P{2});
        % T{2}=P{2};
        % for i=1:8/2-1
        %     T{2}=blkdiag(T{2},P{2});
        % end
        % 
        % T{3}=blkdiag(P{3},P{3});
        % T{4}=P{4};

for k=1:4
%     figure
for i=1:8
%     z{i}=T{1}(i,:)*(y{i}');
    a{i}=T{k}(i,:);
    b{i}=y{i};
    c{i}=a{i}(:)*b{i};
end
for j=1:8
for i=1:8
%     hold on
%     subplot(8,1,j),plot(x{i},c{j}(i,:),'k','linewidth',2),xticks(linspace(0,1,9))
end
end
end