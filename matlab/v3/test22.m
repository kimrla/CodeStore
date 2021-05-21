clear all

NumP = 128;
dh = 1 / NumP;
t = 2*pi*[dh/2:dh:1]';
PP = (0.95+0.1*rand(NumP,1)).*exp(i*t);
P = [real(PP),imag(PP)];  
k = 3;
N = 4;
% CList = [0,1,0];
CList = [0,1,1;
         0.5,0.5,0;
         0.25,0.25,0;
         0.75,0.75,1;
         0.125,0.125,1;
         0.375,0.375,1;
         0.625,0.625,1;
         0.875,0.875,1];
     
Lambda = LSCurFit_V(P,k,N,CList);

%% Show2
figure,
plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
VCompose(Lambda,k,N)
axis equal
axis off