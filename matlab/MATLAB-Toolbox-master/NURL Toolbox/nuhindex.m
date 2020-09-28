function varargout = nuhindex(m, n, k)

% nuhindex: index of nuh basis
%
% Calling Sequences:
%
%     [pt, jc] = nuhindex (m)
%     [pt, jc, hs] = nuhindex (m, n)
%     [pt, jc, hs, td] = nuhindex (m, n, k)
%
% INPUTS:
%
%      m, n, k - the number of interpolation points on two directions
%
% OUTPUT:
%
%   pt  - index for function values
%   jc  - index for first derivatives (Jacobian).
%   hs - index for second derivatives (Hessian).
%   td - index for third derivatives.
%
% See also:
%  
%     nuhgintvdeval
%

if nargin==1
    pt=false(2*m, 1); jc=pt; 
    pt(1:2:2*m-1)=true; 
    jc(2:2:2*m)=true; 
    varargout{1}=pt; 
    varargout{2}=jc; 
elseif nargin==2
    [pt1, jc1]=nuhindex(m);
    [pt2, jc2]=nuhindex(n);
    pt=logical(kron(pt2, pt1));
    jc{1}=logical(kron(pt2, jc1));
    jc{2}=logical(kron(jc2, pt1));
    hs=logical(kron(jc2, jc1));
    varargout{1}=pt; 
    varargout{2}=jc; 
    varargout{3}=hs; 
elseif nargin==3
    [pt1, jc1]=nuhindex(m);
    [pt2, jc2]=nuhindex(n);
    [pt3, jc3]=nuhindex(k);
    pt=logical(kron(pt3, kron(pt2, pt1)));
    jc{1}=logical(kron(pt3, kron(pt2, jc1)));
    jc{2}=logical(kron(pt3, kron(jc2, pt1)));
    jc{3}=logical(kron(jc3, kron(pt2, pt1)));
    hs{1}=logical(kron(pt3, kron(jc2, jc1)));
    hs{2}=logical(kron(jc3, kron(pt2, jc1)));
    hs{3}=logical(kron(jc3, kron(jc2, pt1)));
    td=logical(kron(jc3, kron(jc2, jc1)));
    varargout{1}=pt; 
    varargout{2}=jc; 
    varargout{3}=hs; 
    varargout{4}=td; 
end





