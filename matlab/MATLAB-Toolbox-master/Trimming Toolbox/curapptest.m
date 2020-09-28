
 coefs =[ 6.0  0.0  6.0  1;
         -5.5  0.5  5.5  1;
         -5.0  1.0 -5.0  1;
          4.5  1.5 -4.5  1;
          4.0  2.0  4.0  1;
         -3.5  2.5  3.5  1;
         -3.0  3.0 -3.0  1;
          2.5  3.5 -2.5  1;
          2.0  4.0  2.0  1;
         -1.5  4.5  1.5  1;
         -1.0  5.0 -1.0  1;
          0.5  5.5 -0.5  1;
          0.0  6.0  0.0  1]';
knots = [0 0 0 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 1 1 1];
crv = nrbmak(coefs,knots);
% crv=nrbtestcrv;
u=linspace(0,1,50);N=50-2;
u=u(2:end-1);
P=nrbeval(crv,u);
delta=(u(2)-u(1))/5;
u1=u+ones(1,N)*delta;
u2=u-ones(1,N)*delta;
points1=nrbeval(crv,u1);  
points2=nrbeval(crv,u2);
points=points1-points2;
mag=sqrt(sum(points.*points));
kp=2*delta./mag;
cp=zeros(3,N);
cp(3,:)=kp*5;
figure;hold on;
nrbplot(crv,150);
quiver3(P(1,:),P(2,:),P(3,:),cp(1,:),cp(2,:),cp(3,:),'LineWidth',2);


