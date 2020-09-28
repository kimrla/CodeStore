%% Check the parameters of the nurbs node curve
% The nurbs node information is as follows:
%   nurbs.nDegree ------ The degree of nurbs node
%   nurbs.vecKnots ----- The knot vector of nurbs node
%   nurbs.vecPoles ----- The poles of nurbs node
%   nurbs.vecWeights --- The weights of nurbs node
%   nurbs.bRational ---- Whether it is a rational B-spline node
function nurbsNew = CheckNurbs(nurbs)
    % Check if it is a nurbs node
    if nurbs.nDegree <= 0 || length(nurbs.vecKnots) < 1 || size(nurbs.vecPoles,1) < 1 ...
            || length(nurbs.vecKnots) ~= size(nurbs.vecPoles,1) + nurbs.nDegree + 1
        error('Wrong parameters: the degree of nurbs node is %d, the number of knot vector is %d, and the number of poles is %d\n',...
            nurbs.nDegree, length(nurbs.vecKnots), size(nurbs.vecPoles,1));
    end
    if length(nurbs.vecWeights) > 1 && length(nurbs.vecWeights) ~= size(nurbs.vecPoles,1)
        error('Wrong parameters: the number of poles is %d, and the number of weights is %d\n',...
            size(nurbs.vecPoles,1), length(nurbs.vecWeights));
    end
    % Check if it is a rational B-spline node
    global g_nCompareError;
    nurbsNew = nurbs;
    nurbsNew.bRational = false;
    for i = 2:length(nurbsNew.vecWeights)
        if abs(nurbsNew.vecWeights(i) - nurbsNew.vecWeights(1)) > g_nCompareError
            nurbsNew.bRational = true;
            break;
        end
    end
end