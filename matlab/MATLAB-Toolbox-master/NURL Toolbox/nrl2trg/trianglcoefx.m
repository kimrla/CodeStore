function ctx=trianglcoefx(S,T)

% Derivatives of coefficients of blending function on triangles

[n,m]=size(S); 
c1=zeros(n,m);
c2=c1; c3=c1; 

% Case S==1
pp=S==1;
c1(pp)=1./(1-T(pp));
c2(pp)=0;
c3(pp)=0;

% Case T==1
qq=T==1;
c1(qq)=0;
c2(qq)=-1;
c3(qq)=T(qq)./(1-S(qq)).^2;

% Otherwise
rr=~(qq | pp);
c1(rr)=1./(1-T(rr));
c2(rr)=0;
c3(rr)=T(rr)./(1-S(rr)).^2;

ctx={c1, c2, c3};






