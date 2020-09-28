function [ic, ik] = nuldegelev(d, c, knt, v, t)

% NULDEGELEV:  Degree elevate of NUL. 
% 
% Calling Sequence: 
% 
%   [ic, ik] = nuldegelev(d, c, k, v, t) 
% 
%   INPUT: 
% 
%   d - Degree of the NUL. 
%   c - Control points, matrix of size (dim,nc). 
%   k - Knot sequence, row vector of size nk. 
%   v - intervals 
%   t - Raise the nuls degree t times. 
% 
%   OUTPUT: 
% 
%   ic - Control points of the new NUL. 
%   ik - Knot vector of the new NUL. 
% 

if t<=0
    return;
end

if abs(t-fix(t))>0
    error('The degree elevated should be an interger');
end

[mc, nc] = size(c); 
mv=length(v); 
ic=zeros(mc, nc+(mv-1)*t); 
ik=zeros(1, nc+(mv-1)*t); 
s=0; 
for i=1:mv-1 
    a=v(i); b=v(i+1); 
    p=knt>=a & knt<=b; 
    k=length(find(p)); 
    kk=d+1+t; 
    if k<kk 
        tt=linspace(a, b, kk); 
        pntsi=nulintvdeval(tt, knt(p), c(:,p), d, 0); 
        q=s+1:s+kk; 
        ic(:,q)=pntsi(:,:,1); 
        ik(q)=tt(:); 
    else
        kk=k; q=s+1 : s+kk; 
        ic(:,q)=c(:,p); 
        ik(q)=knt(p); 
    end    
    s=s+kk-1; 
end

ic(:, s+2:nc+(mv-1)*t)=[]; 
ik(s+2:nc+(mv-1)*t)=[]; 



