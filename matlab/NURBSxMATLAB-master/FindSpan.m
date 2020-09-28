%% Determine the knot span index (left)
% Input:  u, degree, knotVector
% Output: the knot span index
% need to be polished if double knots or more
% Xu Yi, 2019

%%
function knotspanIndex = FindSpan(u, knotVector)
knotNum = length(knotVector);
if u == knotVector(knotNum) % special case
    knotspanIndex = knotNum;
    while u == knotVector(knotspanIndex-1)
        knotspanIndex = knotspanIndex-1;
    end
elseif u == knotVector(1) % special case
    knotspanIndex = 1;
    while u == knotVector(knotspanIndex+1)
        knotspanIndex = knotspanIndex+1;
    end
else % do binary search
    temp_low = 1;
    temp_high = knotNum;
    temp_mid = floor( (temp_low+temp_high) /2);
    while u < knotVector(temp_mid) || u >= knotVector(temp_mid+1)
        if u == knotVector(temp_mid)
            break
        elseif u < knotVector(temp_mid)
            temp_high = temp_mid;
        else
            temp_low = temp_mid;
        end
        temp_mid = floor( (temp_low+temp_high) /2);
    end
    knotspanIndex = temp_mid;
end
end
