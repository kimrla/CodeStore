function W=nuhsrftransbk(W, bv, bt, Qx, Qn)

% uhsrftransbk: transform boundary DOFs back to natrual coodinates
%
% Calling Sequences:
%
%     W=uhsrftransbk(W, bv, bt, Qx)
%
%     W=uhsrftransbk(W, bv, bt, Qx, Qn)
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
%  Discriptions:
%      
%      This rutine is used for vibration anslysis of thin plate
%

if nargin==4
    W=W(bv);
    for i=1:length(bt)
        p=4*(bt(i)-1)+[2, 3];
        W(p)=Qx{i}*W(p);
    end
    W(bv)=W;
elseif nargin==5
    W=W(bv);
    for i=1:length(bt)
        p=4*(bt(i)-1)+[2, 3];
        W(p)=Qx{i}*Qn{i}*W(p);
    end
    W(bv)=W;
end





