%% ========================================================================
function [b, c] = otherDim(a)
% return set [1 2 3] without k
b = mod(a+1,3)+1;  % b and c are vertices which are on the same side of the plane
c = 6-a-b;         % a+b+c = 6
end