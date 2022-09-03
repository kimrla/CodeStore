clear all
x=0:0.001:1;
k=3;
N=4;
V=LSMatrix_V(k,N,x');
V=V';
for i=1:(k+1)*2
    subplot((k+1)*2,1,i),plot(x,V(i,:),'color','k','LineWidth',1.5)
end
figure
for i=(k+1)*2+1:(k+1)*2^2
    subplot((k+1)*2,1,i-(k+1)*2),plot(x,V(i,:),'color','k','LineWidth',1.5)
end
figure
for i=1+(k+1)*2^2:(k+1)*2^2+(k+1)*2
    subplot((k+1)*2,1,i-(k+1)*2^2),plot(x,V(i,:),'color','k','LineWidth',1.5)
end