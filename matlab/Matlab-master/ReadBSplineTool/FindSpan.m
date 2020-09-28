%% Query the subscript of the knot in the knot interval
function nIndex = FindSpan(nDegree, vecKnots, nKnot)
    % Check parameters
    global g_nCompareError;
    nLength = length(vecKnots);
    if nKnot < vecKnots(nDegree + 1) || nKnot > vecKnots(nLength - nDegree)
        if abs(nKnot - vecKnots(nDegree + 1)) < g_nCompareError
            nKnot = vecKnots(nDegree + 1);
        elseif abs(nKnot - vecKnots(nLength - nDegree)) < g_nCompareError
            nKnot = vecKnots(nLength - nDegree);
        else
            error('Knot value u=%f exceeds the range of knot vector: lower bound %f, upper bound %f\n',...
                nKnot, vecKnots(nDegree + 1), vecKnots(nLength - nDegree));
        end
    end
    % Special case
    if nKnot == vecKnots(nLength - nDegree)
        nIndex = nLength - nDegree - 1;
        while vecKnots(nIndex) == vecKnots(nLength - nDegree)
            nIndex = nIndex - 1;
        end
    else
        % Binary search
        nLow = nDegree + 1;
        nHigh = nLength - nDegree;
        nMid = floor((nLow + nHigh) / 2);
        while nKnot < vecKnots(nMid) || nKnot >= vecKnots(nMid + 1)
            if nKnot < vecKnots(nMid)
                nHigh = nMid;
            else
                nLow = nMid;
            end
            nMid = floor((nLow + nHigh) / 2);
        end
        nIndex = nMid;
    end
end