function [ic, ik] = nulkntins(d, cpts, knt, intervals, ut)

% NULKNTINS:  Insert knots into a NUL
%
% Calling Sequence:
% 
%   [ic, ik] = nulkntins(d, cpts, knt, intervals, intins)
% 
%  INPUT:
%  
%    d  -  spline degree                 integer
%    cpts  -  control points            double  matrix(mc, nc)      
%    knt  -  knot sequence             double  vector(nc) 
%    intervals  -  intervals               double  vector
%    ut  -  new knots                      double  vector(nu)              
% 
%  OUTPUT:
% 
%    ic - new control points double  matrix(mc, nc + sum(ut) ) 
%    ik - new knot sequence  double  vector(nc + sum(ut))
% 

[mc, nc] = size(cpts); 
nu = numel(ut); 

% Check the knots to be inserted
m=length(intervals); 
if nu~=m-1
    error('The length of the vector of numbers of inserted knots should equal to the number of intervals.');
end
pp=ut<=0;
ut(pp)=0;
pp= abs(ut-fix(ut))>0; 
ut(pp)=round(ut(pp)); 

% Inset knots
s=0;
ic=zeros(mc, nc+sum(ut)); 
ik=zeros(1, nc+sum(ut));
for i=1:m-1 
    a=intervals(i); b=intervals(i+1);
    p=knt>=a & knt<=b;
    k=length(find(p)); 
    if ut(i)>0
        kk=k+ut(i); 
        tt=linspace(a, b, kk); 
        pntsi=nulintvdeval(tt, knt(p), cpts(:,p), d, 0); 
        q=s+1:s+kk;
        ic(:,q)=pntsi(:,:); 
        ik(q)=tt(:); 
    else
        kk=k; q=s+1:s+kk; 
        ic(:,q)=cpts(:,p,1); 
        ik(q)=knt(p); 
    end    
    s=s+kk-1; 
end








