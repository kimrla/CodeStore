function unk_x = unk(n,k,x)
%U 此处显示有关此函数的摘要
%   此处显示详细说明
    if n==0
        unk_x=1;
    elseif n==1
        unk_x=3^0.5*(1-2*x);
    elseif n==2 & k==1
        unk_x=3^0.5*(1-4*x)*(x>=0&x<0.5)+3^0.5*(4*x-3)*(x>=0.5&x<=1);
    elseif n==2 & k==2
        unk_x=(1-6*x)*(x>=0&x<0.5)+(5-6*x)*(x>=0.5&x<=1);
    elseif mod(k,2)~=0
        unk_x=unk(n-1,(k+1)/2,2*x)*(x>=0&x<0.5)+unk(n-1,(k+1)/2,2-2*x)*(x>=0.5&x<=1);
    else
        unk_x=unk(n-1,k/2,2*x)*(x>=0&x<0.5)-unk(n-1,k/2,2-2*x)*(x>=0.5&x<=1);
end

