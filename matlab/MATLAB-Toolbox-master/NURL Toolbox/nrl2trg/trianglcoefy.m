function cty=trianglcoefy(S,T)

% Derivatives of coefficients of blending function on triangles
[n,m]=size(S); c1=zeros(n,m);
c2=c1; c3=c1; 

% Case S==1
pp=S==1;
c1(pp)=0;
c2(pp)=-1;
c3(pp)=S(pp)./(1-T(pp)).^2;

% Case T==1
qq=T==1;
c1(qq)=1./(1-S(qq));
c2(qq)=0;
c3(qq)=0;

% Otherwise
rr=~(qq | pp);
c1(rr)=1./(1-S(rr));
c2(rr)=0;
c3(rr)=S(rr)./(1-T(rr)).^2;

cty={c1, c2, c3};





