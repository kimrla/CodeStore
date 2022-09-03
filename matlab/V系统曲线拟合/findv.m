function [group,i,j] = findv(n,k)
%FINDV 此处显示有关此函数的摘要
%   此处显示详细说明
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
end

