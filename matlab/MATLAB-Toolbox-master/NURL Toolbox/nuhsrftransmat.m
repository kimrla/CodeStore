function [Qx, Qn, bv, bt]=nuhsrftransmat(srf, tbp, nbp)

% nuhsrftransmat: transformation matrices for boundary DOFs of a nurl surface
%
% Calling Sequences:
%
%     [Qx, Qn, bv, bt]=nuhsrftransmat(srf, tbp, nbp)
%
% INPUTS:
%
%      srf - a nurl surface
%
%     tbp - boundary tangent vectors
%
%     nbp  - boundary normal vectors
%
% OUTPUT:
%
%   Qx  - transform from derivatives to Cartesian coodinates
%   
%   Qn - transform from derivatives to normal and tangent directions
%
%   bv  - transform DOFs sequences
% 
%   bt  - index of boundary DOFs
%
%  Discriptions:
%      
%      This rutine is used for vibration anslysis of thin plate
%

m=srf.number(1); n=srf.number(2); 
TN=m*n; TB=2*(m+n)-4;
[~, jac] = nrldeval (srf, srf.knots);
Jacb=squeeze(jac{1}(1,:,:).*jac{2}(2,:,:)-jac{1}(2,:,:).*jac{2}(1,:,:));
t=1;
bv=zeros(1, 4*TN);    
for i=1:m*n
    bv((t-1)*4+1:t*4)=[t, TN+t, 2*TN+t, 3*TN+t];
    t=t+1; 
end
Qx=cell(1, TB); Qn=Qx;
bt=zeros(1, TB); 
t=1; jj=[1, n]; 
for i=1:m
    for k=1:2
        j=jj(k); p=(j-1)*m+i; bt(t)=p;
        qs2x(1, 1)=jac{2}(2,i,j)/Jacb(i,j);
        qs2x(1, 2)=-jac{1}(2,i,j)/Jacb(i,j);
        qs2x(2, 1)=-jac{2}(1,i,j)/Jacb(i,j);
        qs2x(2, 2)=jac{1}(1,i,j)/Jacb(i,j);
        qx2n(1,1)=tbp{1}(1,i); 
        qx2n(1,2)=tbp{1}(2,i); 
        qx2n(2,1)=nbp{1}(1,i); 
        qx2n(2,2)=nbp{1}(2,i); 
        Qx{t}=inv(qs2x);
        Qn{t}=inv(qx2n);
        t=t+1;
    end
end
ii=[1, m];
for k=1:2
    for j=2:n-1
        i=ii(k); p=(j-1)*m+i; bt(t)=p;
        qs2x(1, 1)=jac{2}(2,i,j)/Jacb(i,j);
        qs2x(1, 2)=-jac{1}(2,i,j)/Jacb(i,j);
        qs2x(2, 1)=-jac{2}(1,i,j)/Jacb(i,j);
        qs2x(2, 2)=jac{1}(1,i,j)/Jacb(i,j);
        qx2n(1,1)=tbp{1}(1,i); 
        qx2n(1,2)=tbp{1}(2,i); 
        qx2n(2,1)=nbp{1}(1,i); 
        qx2n(2,2)=nbp{1}(2,i); 
        Qx{t}=inv(qs2x);
        Qn{t}=inv(qx2n);
        t=t+1;
    end
end
[bt, I]=sort(bt);
Qx=Qx(I);
Qn=Qn(I);







