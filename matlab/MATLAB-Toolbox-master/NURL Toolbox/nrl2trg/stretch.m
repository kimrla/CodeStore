% stretch 2D grids to 1D vectors and put in one array or the reverse

function [xx, XX]=stretch(X,m,n)

if nargin==1
    [m,n]=size(X); 
    xx=zeros(m*n,1); 
    XX=zeros(m*n);
    for j=1:n
        for i=1:m
            p=(j-1)*m+i;
            xx(p)=X(j,i);
            XX(p,p)=xx(p);
        end
    end
elseif nargin==2
    s=X; t=m;
    m=length(s); n=length(t);
    S=zeros(m,n); T=S;
    for i=1:m
        for j=1:n
            S(i,j)=s(i); T(i,j)=t(j); 
        end
    end   
    xx=S; XX=T;
elseif nargin==3
    xx=zeros(m,n); 
    XX=zeros(m*n);
    for j=1:n
        for i=1:m
            p=(j-1)*m+i;
            xx(i,j)=X(p); 
            XX(p,p)=X(p);
        end
    end
end






