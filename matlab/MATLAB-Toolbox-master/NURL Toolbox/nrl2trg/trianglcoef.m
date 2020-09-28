function ct=trianglcoef(S,T)

% Coefficients of blending function on triangles

[n,m]=size(S); c1=zeros(n,m);
c2=c1; c3=c1; c4=c1; c5=c1;

% Case S==1
pp=S==1; 
c1(pp)=1;
c2(pp)=S(pp)./(1-T(pp));
c3(pp)=(1-S(pp)-T(pp))./(1-T(pp));
c4(pp)=-(1-S(pp)-T(pp));
c5(pp)=-1;

% Case T==1
qq=T==1;
c1(qq)=(1-S(qq)-T(qq))./(1-S(qq));
c2(qq)=0;
c3(qq)=1;
c4(qq)=-(1-S(qq)-T(qq));
c5(qq)=-S(qq).*(1-S(qq)-T(qq))./(1-S(qq));

% Otherwise
rr=~(qq | pp);
c1(rr)=(1-S(rr)-T(rr))./(1-S(rr));
c2(rr)=S(rr)./(1-T(rr));
c3(rr)=(1-S(rr)-T(rr))./(1-T(rr));
c4(rr)=-(1-S(rr)-T(rr));
c5(rr)=-S(rr).*(1-S(rr)-T(rr))./(1-S(rr));

ct={c1, c2, c3, c4, c5};






