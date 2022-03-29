t=0:0.01:1;
y=(t*10-5).^3+(t*10-5).^2-100*t+5;
scatter(t,y);hold on
k=3;
N=1;
% P=[t;y]';
P=y';
% t=canshuhua(P);
Lambda=LSMatrix_V(k,N,t)\P;
VCompose(Lambda,k,N)