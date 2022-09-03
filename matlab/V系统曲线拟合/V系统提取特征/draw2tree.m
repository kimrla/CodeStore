function draw2tree(k,N,lambda)
%% ��ϵ��lambda����״ͼ

s = 1:k+1;
t = k+2:(k+1)*2;
if N>2
    s=[s,kron(k+2:length(lambda)/2,ones(1,2))];
    t=[t,(k+1)*2+1:length(lambda)];
end

% normlambda(1:k+1)=abs(lambda(1:k+1))/sum(abs(lambda(1:k+1)));
% normlambda(k+2:(k+1)*2)=abs(lambda(k+2:(k+1)*2))/sum(abs(lambda(k+2:(k+1)*2)));
normlambda(1:(k+1)*2)=1;
intervalsize(1:(k+1)*2)=1;
for n=3:N %��ÿһ���������ϵ���ֱ�ȡ����ֵ���һ��
    for j=1:k+1
    normlambda((k+j)*2^(n-2)+1:(k+1+j)*2^(n-2))= ...,
        abs(lambda((k+j)*2^(n-2)+1:(k+1+j)*2^(n-2))) ...,
        /max(abs(lambda((k+j)*2^(n-2)+1:(k+1+j)*2^(n-2))));
    end
    intervalsize((k+1)*2^(n-2)+1:(k+1)*2^(n-1))=1/2^(n-2);
end

G = digraph(s,t);
figure
tree=plot(G,'MarkerSize',4,'Layout','Layered');
set(tree,'NodeCData',normlambda)
set(tree,'MarkerSize',intervalsize*80)
tree.Marker='s';
colormap()
colorbar
end
