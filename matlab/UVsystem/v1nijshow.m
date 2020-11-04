clear all
x=0:0.001:1;
nmax=4;
num=0;
for n=1:nmax    
% n=3;
    for i=1:2 
        if n<3
            for t=1:length(x)
                v1nij_x(t)=v1nij(n,i,1,x(t));
            end
            num=num+1;
            subplot(ceil((2^nmax)^0.5),ceil((2^nmax)^0.5),num),plot(x,v1nij_x);
        else
            for j=1:2^(n-2)
                for t=1:length(x)
                    v1nij_x(t)=v1nij(n,i,j,x(t));
                end
                num=num+1;
                subplot(ceil((2^nmax)^0.5),ceil((2^nmax)^0.5),num),plot(x,v1nij_x);
            end
        end
    end
end
