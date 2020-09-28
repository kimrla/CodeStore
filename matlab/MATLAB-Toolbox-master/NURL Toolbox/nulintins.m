function [ic, ik, iv] = nulintins(d, cpts, knt, intervals, intins)

% NULINTINS:  Insert intervals into a NUL
%
% Calling Sequence:
% 
%   [ic, ik, iv] = nulintins(d, cpts, knt, intervals, intins)
% 
%  INPUT:
%  
%    d  -  spline degree                 integer
%    cpts  -  control points            double  matrix(mc, nc)      
%    knt  -  knot sequence             double  vector(nc) 
%    intervals  -  intervals               double  vector
%    intins  -  new intervals            double  vector(nu)              
% 
%  OUTPUT:
% 
%    ic - new control points double  matrix(mc, nc+ ) 
%    ik - new knot sequence  double  vector(nc + )
%    iv - new interval sequence  double  vector(mv)
% 

[mc, nc] = size(cpts); 
intins  = sort(intins); 
nu = numel(intins); 

% Find the intervals that can be inserted
p=intins>0 & intins<1; 
intins=intins(p); 
if isempty(intins)
    ic=cpts;
    ik=knt;
    iv=intervals;
    return;
end
n=length(intins); 
dm = DistanceMatrix(intins', intervals'); 
pin=zeros(1, n); 
for i=1:n 
    if isempty(find(dm(i,:)<1e-6, 1)) 
        pin(i)=1; 
    end
end
pin=logical(pin); 
intins=intins(pin);  
u=FindUnique(intins');
intins=intins(u);
if isempty(intins)
    ic=cpts;
    ik=knt;
    iv=intervals;
    return;
end

% Insert intervals
s=0; 
ic=zeros(mc, nc+nu*d); 
ik=zeros(1, nc+nu*d);
iv=sort([intervals, intins]);
mv=length(iv); 
for i=1:mv-1 
    a=iv(i); b=iv(i+1);
    p=knt>=a & knt<=b;
    k=length(find(p));
    kk=d+1;
    if k<kk
        tt=linspace(a, b, kk); 
        j=find(intervals<=a, 1, 'last');
        aj=intervals(j); bj=intervals(j+1);
        p=knt>=aj & knt<=bj;
        pntsi=nulintvdeval(tt, knt(p), cpts(:,p), d, 0);
        q=s+1:s+kk;
        ic(:,q)=pntsi(:,:);
        ik(q)=tt(:);
    else
        kk=k; tt=linspace(a, b, kk); 
        pntsi=nulintvdeval(tt, knt(p), cpts(:,p), d, 0);
        q=s+1:s+kk;
        ic(:,q)=pntsi(:,:);
        ik(q)=tt(:);
    end
    s=s+kk-1;
end

ic(:, s+2:nc+nu*d)=[];
ik(s+2:nc+nu*d)=[];







