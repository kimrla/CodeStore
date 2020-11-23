clear all
x=0:0.001:1;
nmax=3;
num=0;
k=1;
for n=1:nmax    
% n=3;
    for i=1:k+1 
        if n<3
            for t=1:length(x)
                vknij_x(t)=vknij(k,n,i,1,x(t));
            end
            num=num+1;
%             subplot(ceil(((k+1)/2*2^nmax)^0.5),ceil(((k+1)/2*2^nmax)^0.5),num),plot(x,vknij_x);
            subplot(1/2*2^nmax,k+1,num),plot(x,vknij_x);
        else
            for j=1:2^(n-2)
                for t=1:length(x)
                    vknij_x(t)=vknij(k,n,i,j,x(t));
                end
                num=num+1;
%                 subplot(ceil(((k+1)/2*2^nmax)^0.5),ceil(((k+1)/2*2^nmax)^0.5),num),plot(x,vknij_x);
            subplot(1/2*2^nmax,k+1,num),plot(x,vknij_x);
            end
        end
        ylabel(['i=',num2str(i)])
    end
    xlabel(['n=',num2str(n)])
end
