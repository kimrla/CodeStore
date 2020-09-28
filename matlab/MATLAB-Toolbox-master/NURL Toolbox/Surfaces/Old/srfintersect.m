function [x, pnts1, pnts2]=srfintersect(srf1, srf2, x)

F=zeros(4,1); dF=zeros(4);
for i=1:5
    [pnts1, jac1, hess1]=nrldeval(srf1, x(1:2)'); 
    [pnts2, jac2, hess2]=nrldeval(srf2, x(3:4)');  
    dr=pnts1-pnts2; 
    F(1)=dot(dr, jac1{1}); 
    F(2)=dot(dr, jac1{2}); 
    F(3)=-dot(dr, jac2{1}); 
    F(4)=-dot(dr, jac2{2}); 
    dF(1,1)=dot(dr, hess1{1,1})+dot(jac1{1}, jac1{1}); 
    dF(1,2)=dot(dr, hess1{1,2})+dot(jac1{1}, jac1{2}); 
    dF(1,3)=-dot(jac1{1}, jac2{1}); 
    dF(1,4)=-dot(jac1{1}, jac2{2}); 
    dF(2,1)=dF(1,2); dF(3,1)=dF(1,3); dF(4,1)=dF(1,4);  
    dF(2,2)=dot(dr, hess1{2,2})+dot(jac1{2}, jac1{2}); 
    dF(2,3)=-dot(jac1{2}, jac2{1}); 
    dF(2,4)=-dot(jac1{2}, jac2{2}); 
    dF(3,2)=dF(2,3); dF(4,2)=dF(2,4); 
    dF(3,3)=-dot(dr, hess2{1,1})+dot(jac2{1}, jac2{1}); 
    dF(3,4)=-dot(dr, hess2{1,2})+dot(jac2{1}, jac2{2}); 
    dF(4,4)=-dot(dr, hess2{2,2})+dot(jac2{2}, jac2{2}); 
    dF(4,3)=dF(3,4); 
    x=x-dF\F; 
    pp=x>1; x(pp)=1; 
    pp=x<0; x(pp)=0; 
end





