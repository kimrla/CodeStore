%% Divide nurbs node into nurbs node segments without heavy knots
% The nurbs node information is as follows:
%   nurbs.nDegree ------ The degree of nurbs node
%   nurbs.vecKnots ----- The knot vector of nurbs node
%   nurbs.vecPoles ----- The poles of nurbs node
%   nurbs.vecWeights --- The weights of nurbs node
%   nurbs.bRational ---- Whether it is a rational B-spline node
function vecNurbs = DivideNurbs(nurbs)
    % Check parameters
    nurbs = CheckNurbs(nurbs);
    % Insert the knot at the inner heavy knot so that the repeat degree of the heavy knot is nurbs.nDegree + 1,
    % and then divide the nurbs node into nurbs node segments without heavy knots
    global g_nCompareError;
    nStartIndex = 1;
    while true
        if abs(nurbs.vecKnots(nStartIndex) - nurbs.vecKnots(nurbs.nDegree+1)) < g_nCompareError
            break;
        end
        nStartIndex = nStartIndex + 1;
    end
    nEndIndex = length(nurbs.vecKnots);
    while true
        if abs(nurbs.vecKnots(nEndIndex) - nurbs.vecKnots(length(nurbs.vecKnots)-nurbs.nDegree)) < g_nCompareError
            break;
        end
        nEndIndex = nEndIndex - 1;
    end
    % Count the knot values that need to be inserted
    vecInsertKnot = nurbs.vecKnots;
    nIndex = 0;
    nFlagKnot = nurbs.vecKnots(nStartIndex);
    nRepeatCount = 1;
    for i = nStartIndex+1:nEndIndex
        if abs(nurbs.vecKnots(i) - nFlagKnot) < g_nCompareError
            nRepeatCount = nRepeatCount + 1;
            continue;
        end
        if nRepeatCount > 1 || nurbs.nDegree == 1
            nInsertCount = nurbs.nDegree + 1 - nRepeatCount;
            for j = 1:nInsertCount
                nIndex = nIndex + 1;
                vecInsertKnot(nIndex) = nFlagKnot;
            end
        end
        nFlagKnot = nurbs.vecKnots(i);
        nRepeatCount = 1;
    end
    if nRepeatCount > 1
        nInsertCount = nurbs.nDegree + 1 - nRepeatCount;
        for j = 1:nInsertCount
            nIndex = nIndex + 1;
            vecInsertKnot(nIndex) = nFlagKnot;
        end
    end
    vecInsertKnot = vecInsertKnot(1:nIndex);
    % Knot refinement
    nurbs = RefineNurbs(nurbs, vecInsertKnot);
    % Divide the nurbs node
    vecNurbs = SectionNurbs(nurbs);
end

%%% Nurbs node: knot refinement
function nurbsNew = RefineNurbs(nurbs, vecInsertKnot)
    if isempty(vecInsertKnot)
        nurbsNew = nurbs;
        return;
    end
    % Non-rational B-spline node
    if ~nurbs.bRational
        nurbsNew = RefineBSpline(nurbs, vecInsertKnot);
        return;
    end
    % Rational B-spline node: transform coordinates into high-dimensional homogeneous coordinates
    vecTempPoles = nurbs.vecPoles;
    [nRow, nColumn] = size(vecTempPoles);
    nurbs.vecPoles = zeros(nRow, nColumn+1);
    for i = 1:nRow
        nurbs.vecPoles(i,1:nColumn) = vecTempPoles(i,:) * nurbs.vecWeights(i);
        nurbs.vecPoles(i,nColumn+1) = nurbs.vecWeights(i);
    end
    nurbsNew = RefineBSpline(nurbs, vecInsertKnot);
    % Transform homogeneous coordinates
    vecTempPoles = nurbsNew.vecPoles;
    [nRow, nColumn] = size(vecTempPoles);
    nurbsNew.vecPoles = zeros(nRow, nColumn-1);
    for i = 1:nRow
        nurbsNew.vecPoles(i,:) = vecTempPoles(i,1:nColumn-1) / vecTempPoles(i,nColumn);
    end
    nurbsNew.vecWeights = vecTempPoles(:,nColumn);
end

%%% B-spline node: knot refinement
function nurbsNew = RefineBSpline(nurbs, vecInsertKnot)
    nLength = length(vecInsertKnot);
    % The degree of B-spline node
    nurbsNew.nDegree = nurbs.nDegree;
    % Knot index
    nMinIndex = FindSpan(nurbs.nDegree, nurbs.vecKnots, vecInsertKnot(1));
    nMaxIndex = FindSpan(nurbs.nDegree, nurbs.vecKnots, vecInsertKnot(nLength)) + 1;
    % Copy invariant knot vector
    nurbsNew.vecKnots = zeros(length(nurbs.vecKnots)+nLength,1);
    nurbsNew.vecKnots(1:nMinIndex) = nurbs.vecKnots(1:nMinIndex);
    nurbsNew.vecKnots(nMaxIndex+nLength:end) = nurbs.vecKnots(nMaxIndex:end);
    % Copy invariant poles
    [nRow, nColumn] = size(nurbs.vecPoles);
    nurbsNew.vecPoles = zeros(nRow+nLength,nColumn);
    nurbsNew.vecPoles(1:nMinIndex-nurbs.nDegree,:) = nurbs.vecPoles(1:nMinIndex-nurbs.nDegree,:);
    nurbsNew.vecPoles(nMaxIndex-1+nLength:end,:) = nurbs.vecPoles(nMaxIndex-1:end,:);
    
    nIndex = nMaxIndex + nurbs.nDegree;
    nFlag = nMaxIndex + nurbs.nDegree + nLength;
    % New poles and knot vector
    for j = nLength:-1:1
        while vecInsertKnot(j) <= nurbs.vecKnots(nIndex) && nIndex > nMinIndex
            nurbsNew.vecPoles(nFlag-nurbs.nDegree-1,:) = nurbs.vecPoles(nIndex-nurbs.nDegree-1,:);
            nurbsNew.vecKnots(nFlag) = nurbs.vecKnots(nIndex);
            nFlag = nFlag - 1;
            nIndex = nIndex - 1;
        end
        nurbsNew.vecPoles(nFlag-nurbs.nDegree-1,:) = nurbsNew.vecPoles(nFlag-nurbs.nDegree,:);
        for k = 1:nurbs.nDegree
            nTemp = nFlag - nurbs.nDegree + k;
            nAlpha = nurbsNew.vecKnots(nFlag+k) - vecInsertKnot(j);
            if nAlpha ~= 0
                nAlpha = nAlpha / (nurbsNew.vecKnots(nFlag+k) - nurbs.vecKnots(nIndex-nurbs.nDegree+k));
                nurbsNew.vecPoles(nTemp-1,:) = nAlpha * nurbsNew.vecPoles(nTemp-1,:) + (1 - nAlpha) * nurbsNew.vecPoles(nTemp,:);
            else
                nurbsNew.vecPoles(nTemp-1,:) = nurbsNew.vecPoles(nTemp,:);
            end
        end
        nurbsNew.vecKnots(nFlag) = vecInsertKnot(j);
        nFlag = nFlag - 1;
    end
    % The weights of poles
    nurbsNew.vecWeights = [];
end

%%% Divide the nurbs node
function vecNurbs = SectionNurbs(nurbs)
    vecNurbs = cell(size(nurbs.vecPoles,1),1);
    nIndex = 0;
    % The first nurbs node segment
    global g_nCompareError;
    nStartIndex = 1;
    while true
        if abs(nurbs.vecKnots(nStartIndex) - nurbs.vecKnots(nurbs.nDegree+1)) < g_nCompareError
            break;
        end
        nStartIndex = nStartIndex + 1;
    end
    nEndIndex = length(nurbs.vecKnots);
    while true
        if abs(nurbs.vecKnots(nEndIndex) - nurbs.vecKnots(length(nurbs.vecKnots)-nurbs.nDegree)) < g_nCompareError
            break;
        end
        nEndIndex = nEndIndex - 1;
    end
    nPointIndex = nStartIndex-1;
    nurbsNew.nDegree = nurbs.nDegree;
    nurbsNew.vecKnots = [];
    nurbsNew.vecWeights = [];
    nRepeatCount = 1;
    nFlagKnot = nurbs.vecKnots(nStartIndex);
    for i = nStartIndex+1:nEndIndex
        % Heavy knot
        if abs(nurbs.vecKnots(i) - nFlagKnot) < g_nCompareError
            nRepeatCount = nRepeatCount + 1;
            continue;
        end
        % Inner knot is not duplicated
        if nRepeatCount == 1
            nurbsNew.vecKnots = [nurbsNew.vecKnots; nFlagKnot];
            nFlagKnot = nurbs.vecKnots(i);
            continue;
        end
        % Inner knot repeats: the repeatability is nurbs.nDegree + 1
        if nRepeatCount < nurbs.nDegree + 1
            error('Wrong parameter: knot repeatability is %d\n', nRepeatCount);
        end
        nurbsNew.vecKnots = [nurbsNew.vecKnots; ones(nRepeatCount,1) * nFlagKnot];
        if length(nurbsNew.vecKnots) == nRepeatCount
            nFlagKnot = nurbs.vecKnots(i);
            nRepeatCount = 1;
            continue;
        end
        nPointCount = length(nurbsNew.vecKnots) - nurbsNew.nDegree - 1;
        nurbsNew.vecPoles = nurbs.vecPoles(nPointIndex+1:nPointIndex+nPointCount,:);
        if ~isempty(nurbs.vecWeights)
            nurbsNew.vecWeights = nurbs.vecWeights(nPointIndex+1:nPointIndex+nPointCount,:);
        end
        nPointIndex = nPointIndex + nPointCount;
        nIndex = nIndex + 1; 
        vecNurbs{nIndex} = nurbsNew;
        % Reset identification value
        nurbsNew.vecKnots = ones(nRepeatCount,1) * nFlagKnot;
        nFlagKnot = nurbs.vecKnots(i);
        nRepeatCount = 1;
    end
    % The last nurbs node segment
    nurbsNew.vecKnots = [nurbsNew.vecKnots; ones(nRepeatCount,1) * nFlagKnot];
    if nPointIndex ~= size(nurbs.vecPoles,1)
        nPointCount = length(nurbsNew.vecKnots) - nurbsNew.nDegree - 1;
        nurbsNew.vecPoles = nurbs.vecPoles(nPointIndex+1:nPointIndex+nPointCount,:);
        if ~isempty(nurbs.vecWeights)
            nurbsNew.vecWeights = nurbs.vecWeights(nPointIndex+1:nPointIndex+nPointCount,:);
        end
        nPointIndex = nPointIndex + nPointCount;
        nIndex = nIndex + 1; 
        vecNurbs{nIndex} = nurbsNew;
    end
    vecNurbs = vecNurbs(1:nIndex);
    if nPointIndex ~= size(nurbs.vecPoles,1) - (length(nurbs.vecKnots) - nEndIndex)
        error('Wrong result: Missing poles!\n');
    end
end