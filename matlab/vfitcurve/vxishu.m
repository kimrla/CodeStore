function [V,lamda] = vxishu(k,N,M,t,d)
%VXISHU �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

for m=1:M
    for n=1:N
        %��N��V�������ǵ�n��ĵ�i��ĵ�j��V�� ����n=log2(N/(k+1))+1����ȡ��, i=N-(k+1)*2^(n-2)��2^(n-2)��������ȡ����ȡ��Ϊj
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

