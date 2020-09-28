function Ke=nuhsrfmcondens(Ke, Qc, tn)

p=1:3*tn; q=3*tn+1:4*tn;
r=2*tn+1:3*tn;
K11=Ke(p, p); K12=Ke(p, q);
K22=Ke(q, q);
Kt=K12(r,:)*Qc;
K11(r,r)=K11(r,r)+Kt+Kt'+Qc'*K22*Qc;
Ke=K11;