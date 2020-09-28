function Ke=nuhsrftrans(Ke, bv, bt, Qx, Qn)

% nuhsrftrans: transform the boundaary DOFs of stiffness or mass matrix
%
% Calling Sequences:
%
%     Ke=nuhsrftrans(Ke, bv, bt, Qx)
%
%     Ke=nuhsrftrans(Ke, bv, bt, Qx, Qn)
%
% INPUTS:
%
%     Ke - stiffness or mass matrix before transformation
%
%     bv  - transform DOFs sequences
% 
%     bt  - index of boundary DOFs
%
%     Qx  - transform from derivatives to Cartesian coodinates
%   
%     Qn - transform from derivatives to normal and tangent directions
%
% OUTPUT:
%
%     Ke - stiffness or mass matrix after transformation
%
%  Discriptions:
%      
%      This rutine is used for vibration anslysis of thin plate
%

if nargin==4
    Ke=Ke(bv, bv);
    for i=1:length(bt)
        p=4*(bt(i)-1)+[2, 3];
        Ke(:, p)=Ke(:, p)*Qx{i};
        Ke(p, :)=Qx{i}'*Ke(p, :);
    end
    Ke(bv, bv)=Ke(:, :); 
elseif nargin==5
    Ke=Ke(bv, bv);
    for i=1:length(bt)
        p=4*(bt(i)-1)+[2, 3];
        Ke(:, p)=Ke(:, p)*Qx{i}*Qn{i};
        Ke(p, :)=Qn{i}'*Qx{i}'*Ke(p, :);
    end
    Ke(bv, bv)=Ke(:, :); 
end







