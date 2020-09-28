%% Calculate the derivative at the specified knot of nurbs node
% The nurbs node information is as follows:
%   nurbs.nDegree ------ The degree of nurbs node
%   nurbs.vecKnots ----- The knot vector of nurbs node
%   nurbs.vecPoles ----- The poles of nurbs node
%   nurbs.vecWeights --- The weights of nurbs node
%   nurbs.bRational ---- Whether it is a rational B-spline node
function [nxDeriv0, nxDeriv1, nxDeriv2] = GetNurbsDeriv(nurbs, nKnot)
    % Check parameters
    nurbs = CheckNurbs(nurbs);
    % Classification processing
    if nurbs.bRational
        % Rational B-spline node: Using Leibniz's formula
        % Zero-order, first-order, and second-order derivatives of the denominator
        [nDenominatorDeriv0, nDenominatorDeriv1, nDenominatorDeriv2] = GetBSplineDeriv(nurbs, nurbs.vecWeights, nKnot);
        % Zero-order, first-order, and second-order derivatives of the numerator
        vecWeightedPoles = nurbs.vecPoles;
        for i = 1:size(nurbs.vecPoles,1)
            vecWeightedPoles(i,:) = nurbs.vecPoles(i,:) * nurbs.vecWeights(i);
        end
        [nxNumeratorDeriv0, nxNumeratorDeriv1, nxNumeratorDeriv2] = GetBSplineDeriv(nurbs, vecWeightedPoles, nKnot);
        % Zero-order derivative (point coordinate)
        nxDeriv0 = nxNumeratorDeriv0 / nDenominatorDeriv0;
        % First-order derivative
        nxDeriv1 = (nxNumeratorDeriv1 * nDenominatorDeriv0 - nxNumeratorDeriv0 * nDenominatorDeriv1)...
            / (nDenominatorDeriv0 * nDenominatorDeriv0);
        % Second-order derivative
        nxDeriv2 = (nxNumeratorDeriv2 * nDenominatorDeriv0 - 2 * nxDeriv1 * nDenominatorDeriv0 * nDenominatorDeriv1...
            - nxNumeratorDeriv0 * nDenominatorDeriv2) / (nDenominatorDeriv0 * nDenominatorDeriv0);
    else
        % Non-rational B-spline node: Solving by De-Boor recursion formula
        [nxDeriv0, nxDeriv1, nxDeriv2] = GetBSplineDeriv(nurbs, nurbs.vecPoles, nKnot);
    end
end

%%% Calculate the derivatives at a specified knot of a non-rational B-spline node
function [nxDeriv0, nxDeriv1, nxDeriv2] = GetBSplineDeriv(nurbs, vecWeightedPoles, nKnot)
    % Get valid poles
    [vecPoles, nKnotIndex] = GetValidPoles(nurbs, vecWeightedPoles, nKnot);
    % Solve using De-Boor recursive algorithm
    vecTempPoles = vecPoles;
    deBoorInfo.nKnot = nKnot;
    deBoorInfo.nKnotIndex = nKnotIndex;
    deBoorInfo.nStartIndex = 1;
    deBoorInfo.nEndIndex = nurbs.nDegree - 2;
    [vecPoles, vecTempPoles] = GetIteratePole(nurbs, deBoorInfo, vecPoles, vecTempPoles);
    if nurbs.nDegree > 1
        % Second-order derivative
        nxDeriv2 = nurbs.nDegree * (nurbs.nDegree - 1) / (nurbs.vecKnots(nKnotIndex+1) - nurbs.vecKnots(nKnotIndex))...
            * ((vecPoles(nurbs.nDegree+1, :) - vecPoles(nurbs.nDegree, :)) / (nurbs.vecKnots(nKnotIndex+2)...
            - nurbs.vecKnots(nKnotIndex)) - (vecPoles(nurbs.nDegree, :) - vecPoles(nurbs.nDegree-1, :))...
            / (nurbs.vecKnots(nKnotIndex+1) - nurbs.vecKnots(nKnotIndex-1)));
    else
        nxDeriv2 = zeros(size(vecPoles(1,:)));
    end
    
    deBoorInfo.nStartIndex = deBoorInfo.nEndIndex + 1;
    deBoorInfo.nEndIndex = deBoorInfo.nStartIndex;
    [vecPoles, vecTempPoles] = GetIteratePole(nurbs, deBoorInfo, vecPoles, vecTempPoles);
    % First-order derivative
    nxDeriv1 = nurbs.nDegree / (nurbs.vecKnots(nKnotIndex+1) - nurbs.vecKnots(nKnotIndex)) * (vecPoles(nurbs.nDegree+1,:)...
        - vecPoles(nurbs.nDegree,:));
    
    deBoorInfo.nStartIndex = deBoorInfo.nEndIndex + 1;
    deBoorInfo.nEndIndex = deBoorInfo.nStartIndex;
    [vecPoles, ~] = GetIteratePole(nurbs, deBoorInfo, vecPoles, vecTempPoles);
    % Zero-order derivative (point coordinate)
    nxDeriv0 = vecPoles(nurbs.nDegree+1, :);
end

%%% Get valid poles
function [vecValidPoles, nKnotIndex] = GetValidPoles(nurbs, vecOrigPoles, nKnot)
    % Query the subscript of the knot in the knot interval
    nKnotIndex = FindSpan(nurbs.nDegree, nurbs.vecKnots, nKnot);
    % Poles for non-zero B-spline basis functions
    vecValidPoles = vecOrigPoles(nKnotIndex-nurbs.nDegree:nKnotIndex,:);
end

%%% Get iterative poles
function [vecPolesNew, vecTempPolesNew] = GetIteratePole(nurbs, deBoorInfo, vecPoles, vecTempPoles)
    vecPolesNew = vecPoles;
    vecTempPolesNew = vecTempPoles;
    % De-Boor recursive algorithm calculates iterative poles
    for i = max(deBoorInfo.nStartIndex,1):deBoorInfo.nEndIndex
        for j = i:nurbs.nDegree
            nTempIndex = deBoorInfo.nKnotIndex - nurbs.nDegree + j;
            nAlpha = (deBoorInfo.nKnot - nurbs.vecKnots(nTempIndex)) / (nurbs.vecKnots(nurbs.nDegree+1+nTempIndex-i)...
                - nurbs.vecKnots(nTempIndex));
            vecTempPolesNew(j+1,:) = (1 - nAlpha) * vecPolesNew(j,:) + nAlpha * vecPolesNew(j+1,:);
        end
        vecPolesNew = vecTempPolesNew;
    end
end