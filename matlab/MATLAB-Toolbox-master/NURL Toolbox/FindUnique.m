function [u, left]=FindUnique(vectmat, tol)

% Find unique elements for vector matrix
% 
% Input : 
% 
%    vectmat - a (n*d) vector matrix where 
%                    n is the number of elements and
%                    d is the dimension of the vector
%    tol - tollerance
%  
% Output :
% 
%    u - indexes of the unique elements
%    left - indexes of the duplicated elements
%             corresponding to u
% 

[n,~]=size(vectmat);
if n==1
    u=1; left=[];
    return;
end

dm = DistanceMatrix(vectmat, vectmat);
if nargin==1
    tol=max(max(dm))*1e-6;
end
uu=ones(1,n); uniq0=cell(n,1);
dm=(dm<tol);
for i=1:n
    uniq0{i}=find(dm(i,:));
end
for i=1:n    
    for j=i+1:n
        if dm(i,j)==1
            dm(:,j)=0;
        end
    end
    if dm(i,i)==0
        uu(i)=0;
    end
end
uu=uu==1;

u=find(uu)'; n=length(u);
left=cell(n,1);
for i=1:n
    tt=uniq0{u(i)}~=u(i);
    left{i}=uniq0{u(i)}(tt);
end



