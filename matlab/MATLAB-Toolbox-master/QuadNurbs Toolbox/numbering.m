function [nums, C, dump, uniq]=numbering(A, tol)
%
% numbering: Get node numbers of a sets of nodes
% 
% Calling Sequences:
% 
%     [nums, C, dump, uniq]=numbering(A)
% 
%     [nums, C, dump, uniq]=numbering(A, tol)
% 
% INPUTS:
%
%      A - A vector or a m*dim array, where dim is the dimension.
%
%      tol - Tollerance. Default value is 1e-6 multiply  
%              the minimum distance of the nodes.
%
% OUTPUT:
% 
%      nums - Global node numbers of each node. A=C(nums, :).
%
%      C  - Unique values in the array A. 
%
%      dump -  Indexes correspond to each unique node. 
%                    dump{k} includes all dumplicated nodes at C(k, :).
%
%      uniq  -   Indexes of unique nodes. C=A(uniq, :).
%
%  Description:
%
%      See also UniqueId, where [C, ia, ic, rep]=UniqueId(A). 
%      Here ia = uniq, ic = nums, rep = dump.
%

if nargin==1
    [uniq, left]=FindUnique(A);
elseif nargin==2
    [uniq, left]=FindUnique(A, tol);
end

nums=zeros(size(A,1), 1);
C=A(uniq,:);
TN=length(uniq);
dump=cell(TN,1);
for i=1:TN
    dump{i}=[uniq(i), left{i}];
    nums([uniq(i), left{i}])=i;
end





