clear; clc;

order=3; 
knots=linspace(-1, 1, 3);

[x, C]=GaussP(knots, order+1);

plot(x,0*x,'ro');  title('Gauss points');
axis([-1,1,-1,1]);

f=sqrt(x+1.5); ef=2.399529;
intf=sum(C.*f);

syms s;
fs=s+2*s.^2+s.^7+s.^6;
f=x+2*x.^2+x.^7+x.^6;
efs=double(int(fs,-1,1));
intfs=sum(C.*f);




