function lg=isderiv(dnrb)
% Determine if the 1st or 2nd derivative exists, meaning there's no NAN or
% INF of the element in the dnrb.coefs, see [dnrb,dnrb2]=nrbderiv(nrb).
% If the input is tangent vector=[0,0,0] or [inf], then lg=1.

% Input:
%   dnrb: 1st or 2nd derivative structure, see nrbderiv.m.

% Output:
%   lg: Logical scalar where 1: The derivative structure contains NAN or
%       INF, which needs to be replaced by difference quotient, and 0: The
%       derivative is normal.

lg=false;
if iscell(dnrb)   
    if any (arrayfun(@(x) any(isnan(x.coefs(:)) | isinf(x.coefs(:))),dnrb))
        lg=true;
        return;
    end
elseif isstruct(dnrb)
    if any (any(isnan(dnrb.coefs(:)) | isinf(dnrb.coefs(:))))
        lg=true;
        return;
    end
elseif numel(dnrb)==3 % input equals [0,0,0] or [inf]
    n=norm(dnrb);
    if n<=1e-10 || isinf(n)
        lg=true;
        return;
    end
end