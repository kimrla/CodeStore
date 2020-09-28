function [S, T, Se, Te]=TrigLobatto(t)
% 
% TrigLobatto: Generate non-equally spaces nodes on a unite right triangle
% 
% Calling Sequence:
% 
%     [S, T]=TrigLobatto(t)
%
%     [S, T, Se, Te]=TrigLobatto(t)
% 
% INPUT:
% 
%     t :  An unequally space node vector defined on [0, 1]
% 
% OUTPUT:
%
%    S, T :  Natrual coordinates of unequally spaced nodes 
%            on a unit right triangle.
% 
%    Se, Te :  Natrual coordinates of equally spaced nodes 
%            on a unit right triangle.
% 

m=length(t);
s=linspace(0, 1, m); 
N=m*(m+1)/2;

% Equal nodes on unit right triangle
if nargout>2
    Se=zeros(N, 1); Te=Se;
end
Is=zeros(N, 1); It=Is;
p=1;
for j=1:m
    for i=1:(m-j)+1
        if nargout>2
            Se(p)=s(i); Te(p)=s(j); 
        end
        Is(p)=i; It(p)=j;
        p=p+1;
    end
end

% Transform to be non-equally spaced
Ir=m+2-Is-It;
Sn=t(Is); Tn=t(It); Rn=t(Ir); 
Pn=9*Sn.*Tn.*Rn;
Dn=Sn+Tn+Rn+Pn;
S=(Sn+Pn/3)./Dn; T=(Tn+Pn/3)./Dn; 






