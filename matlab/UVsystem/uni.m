function uni_x = uni(n,i,x)
%U 此处显示有关此函数的摘要
%   此处显示详细说明
    if n==0
        uni_x=1;
    elseif n==1
        uni_x=3^0.5*(1-2*x);
    elseif n==2 
        if i==1
        uni_x=3^0.5*(1-4*x)*(x>=0&x<0.5)+3^0.5*(4*x-3)*(x>=0.5&x<=1);
        elseif  i==2            
        uni_x=(1-6*x)*(x>=0&x<0.5)+(5-6*x)*(x>=0.5&x<=1);
        end
    elseif mod(i,2)~=0
        uni_x=uni(n-1,(i+1)/2,2*x)*(x>=0&x<0.5)+uni(n-1,(i+1)/2,2-2*x)*(x>=0.5&x<=1);
    else
        uni_x=uni(n-1,i/2,2*x)*(x>=0&x<0.5)-uni(n-1,i/2,2-2*x)*(x>=0.5&x<=1);
    end
end

