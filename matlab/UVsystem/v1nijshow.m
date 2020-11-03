clear all
x=0:0.001:1;
for n=1:10    
    for t=1:length(x)
        v1nij_x(t)=v1nij(n,i,j,x(t))