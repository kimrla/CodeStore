function [C, ia, ic, rep]=UniqueId(A)
%
% UniqueId: Find unique values in array and indexes of data with repetitions
% 
% Calling Sequence:
% 
%   C=UniqueId(A)
% 
%   [C, ia, ic, rep]=UniqueId(A)
% 
% INPUT:
% 
%    A :  A vector or a m*dim array, where dim is the dimension.
% 
% OUTPUT:
%
%   C : Unique values in the array.
%
%   ia - Index to A, returned as a column vector of indices 
%         to the first occurrence of repeated elements.
%         C = A(ia) (or C = A(ia,:).
%
%   ic - Index to C, returned as a column vector. ic contains indices, 
%        such that A = C(ic) (or A = C(ic,:) for the 'rows' option).
%
%   rep - Indexes of each data in the array with repetitions.
%
%  Description:
%
%      See also numbering, where [nums, C, dump, uniq] = numbering(A, tol).
%      Here ia = uniq, ic = nums, rep = dump.
% 

[m, n]=size(A);
if m==1 || n==1
    [C,ia,ic] = unique(A);
    if nargout==4
        rep=cell(length(C), 1);
        for j=1:length(A)
            rep{ic(j)}=[rep{ic(j)}, j];
        end
    end
else
    [C,ia,ic] = unique(A, 'rows');
    if nargout==4
        rep=cell(size(C, 1), 1);
        for j=1:size(A, 1)
            rep{ic(j)}=[rep{ic(j)}, j];
        end
    end
end





