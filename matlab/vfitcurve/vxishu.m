function [V,lamda] = vxishu(k,N,M,t,d)
%VXISHU 此处显示有关此函数的摘要
%   此处显示详细说明

for m=1:M
    for n=1:N
        %第N个V基函数是第n组的第i类的第j个V基 其中n=log2(N/(k+1))+1向上取整, i=N-(k+1)*2^(n-2)除2^(n-2)的商向上取整，取余为j
        if n<=k+1
            group=1;
            i=n;
            j=1;
        else
        group=ceil(log2(n/(k+1))+1);
        i=ceil((n-(k+1)*2^(group-2))/2^(group-2));
        j=mod((n-(k+1)*2^(group-2)),2^(group-2));
        if j==0
            j=2^(group-2);
        end
        end
        V(m,n)=vknij(k,group,i,j,t(m));
        
    end
end
lamda=pinv(transpose(V)*V)*transpose(V)*d;
end

