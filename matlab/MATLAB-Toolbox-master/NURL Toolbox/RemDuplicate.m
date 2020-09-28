function [vectmat, uu]=RemDuplicate(vectmat, tol)

% Remove duplicate elements for a (n*d) vector matrix
%    where n is the number of elements and
%    d is the dimension of the vector
%删除重复的元素，删除后按照原有顺序输出，不会由大到小排列顺序。输入是列向量
if nargin==1
    tol=1e-3;
end

dm = DistanceMatrix(vectmat, vectmat);
[n, ~]=size(vectmat);
uu=ones(1,n);
dm=(dm<tol);
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
vectmat=vectmat(uu,:);

