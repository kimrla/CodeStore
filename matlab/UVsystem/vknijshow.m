clear all
x=0:0.001:1;
nmax=4;
num=0;
k=2;
for n=1:nmax    
% n=3;
    for i=1:k+1 
        if n<3
            for t=1:length(x)
                vknij_x(t)=vknij(k,n,i,1,x(t));
            end
            num=num+1;
            subplot(ceil(((k+1)/2*2^nmax)^0.5),ceil(((k+1)/2*2^nmax)^0.5),num),plot(x,vknij_x);
        else
            for j=1:2^(n-2)
                for t=1:length(x)
                    vknij_x(t)=vknij(k,n,i,j,x(t));
                end
                num=num+1;
                subplot(ceil(((k+1)/2*2^nmax)^0.5),ceil(((k+1)/2*2^nmax)^0.5),num),plot(x,vknij_x);
            end
        end
    end
end
