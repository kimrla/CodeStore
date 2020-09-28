%% Scatter nurbs node by specified error
% The nurbs node information is as follows:
%   nurbs.nDegree ------ The degree of nurbs node
%   nurbs.vecKnots ----- The knot vector of nurbs node
%   nurbs.vecPoles ----- The poles of nurbs node
%   nurbs.vecWeights --- The weights of nurbs node
%   nurbs.bRational ---- Whether it is a rational B-spline node
function vecPoints = ScatterNurbs(nurbs, nDeflection)
    % Numerical comparison accuracy
    global g_nCompareError;
    g_nCompareError = 1e-9;
    if nurbs.nDegree == 1
        vecTempPoints = nurbs.vecPoles;
    else
        % Divide nurbs node
        vecTempPoints = nurbs.vecPoles(1,:);
        vecNurbs = DivideNurbs(nurbs);
        for i = 1:length(vecNurbs)
            % Check parameters
            nurbsNew = CheckNurbs(vecNurbs{i});
            % Collect unique knot values in the knot vector
            vecValidKnots = nurbsNew.vecKnots(nurbsNew.nDegree+1:length(nurbsNew.vecKnots)-nurbsNew.nDegree);
            vecDispKnots = vecValidKnots;
            nIndex = 1;
            for j = 2:length(vecValidKnots)
                if abs(vecDispKnots(nIndex) - vecValidKnots(j)) > g_nCompareError / 2
                    nIndex = nIndex + 1;
                    vecDispKnots(nIndex) = vecValidKnots(j);
                end
            end
            vecDispKnots = vecDispKnots(1:nIndex);
            % Divide each two non-repeating knot values into two parts and calculate the scatter points of the nurbs node
            for j = 2:length(vecDispKnots)
                nMidKnot = 0.5 * (vecDispKnots(j-1) + vecDispKnots(j));
                vecFrontPoints = PerformNurbs(nurbsNew, nDeflection, vecDispKnots(j-1), nMidKnot);
                vecBackPoints = PerformNurbs(nurbsNew, nDeflection, nMidKnot, vecDispKnots(j));
                vecTempPoints = [vecTempPoints(1:end-1,:); vecFrontPoints(1:end-1,:); vecBackPoints];
            end
        end
    end
    % Merge collinear points
    vecPoints = zeros(size(vecTempPoints));
    vecPoints(1,:) = vecTempPoints(1,:);
    nIndex = 1;
    nTempIndex = 2;
    while nTempIndex <= size(vecTempPoints,1)
        if norm(vecTempPoints(nTempIndex,:) - vecPoints(nIndex,:)) > g_nCompareError
            nIndex = nIndex + 1;
            vecPoints(nIndex,:) = vecTempPoints(nTempIndex,:);
            break;
        end
        nTempIndex = nTempIndex + 1;
    end
    if nTempIndex >= size(vecTempPoints,1)
        vecPoints = vecPoints(1:nIndex,:);
        return;
    end
    for i = nTempIndex+1:size(vecTempPoints,1)
        if JudgeCollinearPoints(vecPoints(nIndex-1,:), vecPoints(nIndex,:), vecTempPoints(i,:))
            vecPoints(nIndex,:) = vecTempPoints(i,:);
        else
            nIndex = nIndex + 1;
            vecPoints(nIndex,:) = vecTempPoints(i,:);
        end
    end
    vecPoints = vecPoints(1:nIndex,:);
end

%%% Scatter nurbs node at specified interval with specified deflection
function vecPoints = PerformNurbs(nurbs, nDeflection, nStartKnot, nEndKnot)
    vecPoints = [];
    nPointMin = 2;
    knotInfoStart.nKnot = nStartKnot;
    knotInfoEnd.nKnot = nEndKnot;
    [knotInfoStart.nxDeriv0, knotInfoStart.nxDeriv1, knotInfoStart.nxDeriv2] = GetNurbsDeriv(nurbs, knotInfoStart.nKnot);
    [knotInfoEnd.nxDeriv0, knotInfoEnd.nxDeriv1, knotInfoEnd.nxDeriv2] = GetNurbsDeriv(nurbs, knotInfoEnd.nKnot);
    vecPoints = [vecPoints; knotInfoStart.nxDeriv0];
    % Recursive calculation
    global g_nCompareError;
    nEps = g_nCompareError * g_nCompareError;
    vecPoints = QuasiFleche(nurbs, nDeflection * nDeflection, knotInfoStart, knotInfoEnd, nPointMin, nEps, vecPoints);
end

%%% Recursively calculate scatter nurbs node with specified accuracy
function vecPointsNew = QuasiFleche(nurbs, nDeflection2, knotInfoStart, knotInfoEnd, nPointMin, nEps, vecPoints)
    vecPointsNew = vecPoints;
    nPointSize = size(vecPointsNew,1);
    nDeltaKnot = knotInfoEnd.nKnot - knotInfoStart.nKnot;
    knotInfoNext = knotInfoEnd;
    if nPointMin > 2
        nDeltaKnot = nDeltaKnot / (nPointMin - 1);
        knotInfoNext.nKnot = knotInfoStart.nKnot + nDeltaKnot;
        [knotInfoNext.nxDeriv0, knotInfoNext.nxDeriv1, knotInfoNext.nxDeriv2] = GetNurbsDeriv(nurbs, knotInfoNext.nKnot);
    end
    
    nNorme = (norm(knotInfoNext.nxDeriv0 - knotInfoStart.nxDeriv0))^2;
    nTheFleche = 0;
    bFlecheOk = false;
    if nNorme > nEps
        nNorme1 = (norm(knotInfoStart.nxDeriv1))^2;
        nNorme2 = (norm(knotInfoNext.nxDeriv1))^2;
        if nNorme1 > nEps && nNorme2 > nEps
            nxStartNormalized = knotInfoStart.nxDeriv1 / norm(knotInfoStart.nxDeriv1);
            nxNextNormalized = knotInfoNext.nxDeriv1 / norm(knotInfoNext.nxDeriv1);
            nNormeDiff = (norm(nxStartNormalized - nxNextNormalized))^2;
            if nNormeDiff > nEps
                nTheFleche = nNormeDiff * nNorme / 64;
                bFlecheOk = true;
            end
        end
    end
    
    if ~bFlecheOk
        nxMid = 0.5 * (knotInfoNext.nxDeriv0 + knotInfoStart.nxDeriv0);
        [nxVerif, ~, ~] = GetNurbsDeriv(nurbs, knotInfoStart.nKnot + 0.5 * nDeltaKnot);
        nTheFleche = (norm(nxMid - nxVerif))^2;
    end
    
    if nTheFleche < nDeflection2
        vecPointsNew = [vecPointsNew; knotInfoNext.nxDeriv0];
    else
        vecPointsNew = QuasiFleche(nurbs, nDeflection2, knotInfoStart, knotInfoNext, 3, nEps, vecPointsNew);
    end
    
    if nPointMin > 2
        nPointMin = nPointMin - (size(vecPointsNew,1) - nPointSize);
        vecPointsNew = QuasiFleche(nurbs, nDeflection2, knotInfoNext, knotInfoEnd, nPointMin, nEps, vecPointsNew);
    end
end

%%% Check if the three points are in line
function bCollinear = JudgeCollinearPoints(nxPoint1, nxPoint2, nxPoint3)
    global g_nCompareError;
    bCollinear = abs(norm(nxPoint2 - nxPoint1) + norm(nxPoint3 - nxPoint2) - norm(nxPoint3 - nxPoint1)) < g_nCompareError;
end