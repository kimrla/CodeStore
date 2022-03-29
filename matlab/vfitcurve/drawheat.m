function drawheat(k,N,lambda)
%% 热力图
heatlambda=zeros(length(lambda)/2,N);
normheatlambda=heatlambda;
heatlambda(:,1:2)=kron(reshape(lambda(1:(k+1)*2),[],2),ones(2^(N-2),1));
for n=3:N
    heatlambda(:,n)=kron(lambda((k+1)*2^(n-2)+1:(k+1)*2^(n-1)),ones(2^(N-n),1));
end
heatlambda=abs(heatlambda);
normheatlambda(:,1:2)=1;
for n=3:N  %将每一类基函数的系数分别取绝对值后归一化
    for j=1:k+1
        normheatlambda(1+(j-1)*2^(N-2):j*2^(N-2),n)= ...,
            heatlambda(1+(j-1)*2^(N-2):j*2^(N-2),n)/ ...,
            max(abs(lambda((k+j)*2^(n-2)+1:(k+1+j)*2^(n-2))));
    end
end
figure
for j=1:k+1
ylable=[0 1];
xlable=1:N;
subplot(k+1,1,j),imagesc(xlable,ylable,normheatlambda(1+(j-1)*2^(N-2):j*2^(N-2),:))
set(gca,'xtick',1:N)
end

end

