%% Scatter arc into small line segments with specified accuracy
% Parameter strG17: Coordinate plane
% Parameter strCurGId: G command (G2 or G3)
% Parameter nxLastPos: Start point of arc
% Parameter nxCurPos: End point of arc
% Parameter nxCurParam: Arc parameters
% Return value vecPoints: The scatter points of arc
function vecPoints = ScatterArc(strG17, strCurGId, nxLastPos, nxCurPos, nxCurParam, nDeflection)
    % Arc command G02, G03 (G02/G03 X_Y_Z_R_(I_J_K_))
    % G02 specifies clockwise interpolation, G03 is counterclockwise
    % A negative value of R indicates that the arc segment is larger than a semicircle, while a positive value of R
    % indicates that the arc segment is less than or equal to a semicircle.
    % G17 sets the X-Y work plane, G18 sets the Z-X work plane, and G19 sets the Y-Z work plane
    % Parameter check
    if ~strcmp(strCurGId, 'G2') && ~strcmp(strCurGId, 'G3')
        error('Wrong parameter value: %s, this command is not a arc command!', strCurGId);
    end
    % Default parameters
    if isempty(strG17)
        strG17 = 'G17';
    end
    % Convert coordinates to plane coordinates
    [nxLastTwoAxis, nLastThreeAxis] = GetTwoAxisByPlane(strG17, nxLastPos);
    [nxCurTwoAxis, nCurThreeAxis] = GetTwoAxisByPlane(strG17, nxCurPos);
    if nLastThreeAxis ~= nCurThreeAxis
        error('The start and end points of the arc are not in the same coordinate plane!');
    end
    % Arc center
    [nxArcCenter, ~] = GetTwoAxisByPlane(strG17, nxCurParam(1:3));
    global g_IJLIncrementalMode;
    if g_IJLIncrementalMode
        nxArcCenter = nxArcCenter + nxLastTwoAxis;
    end
    % Unified radius programming to center point programming
    if ~isinf(nxCurParam(4))
        nxArcCenter = TransRadiusToArcCenter(strCurGId, nxLastTwoAxis, nxCurTwoAxis, nxCurParam(4));
    end
    % Calculate the angle
    nxStartVec = (nxLastTwoAxis - nxArcCenter) / norm(nxLastTwoAxis - nxArcCenter);
    nxEndVec = (nxCurTwoAxis - nxArcCenter) / norm(nxCurTwoAxis - nxArcCenter);
    nCosValue = nxStartVec(1) * nxEndVec(1) + nxStartVec(2) * nxEndVec(2);
    nSinValue = nxStartVec(1) * nxEndVec(2) - nxStartVec(2) * nxEndVec(1);
    if strcmp(strCurGId, 'G2')
        nSinValue = -nSinValue;
    end
    if nCosValue > 1
        nCosValue = 1;
    elseif nCosValue < -1
        nCosValue = -1;
    end
    nAngle = acos(nCosValue);
    if nSinValue < 0
        nAngle = 2 * pi - nAngle;
    end
    % Whether it is a full circle
    if 0 == nAngle && isinf(nxCurParam(4))
        nAngle = 2 * pi;
    end
    % Scatter arc based on accuracy
    nRadius = norm(nxLastTwoAxis - nxArcCenter);
    nStepAngle = 2 * acos(1 - min(nDeflection, nRadius) / nRadius);
    nCount = ceil(nAngle / nStepAngle) + 1;
    nStepAngle = nAngle / nCount;
    if strcmp(strCurGId, 'G2')
        nStepAngle = -nStepAngle;
    end
    % Calculate scatter points
    vecPoints = zeros(nCount+1,3);
    nxTempAxis = zeros(1,2);
    for i = 0:nCount
        nStartAngle = i * nStepAngle;
        nxTempAxis(1) = nxArcCenter(1) + nRadius * (nxStartVec(1) * cos(nStartAngle) - nxStartVec(2) * sin(nStartAngle));
        nxTempAxis(2) = nxArcCenter(2) + nRadius * (nxStartVec(1) * sin(nStartAngle) + nxStartVec(2) * cos(nStartAngle));
        vecPoints(i+1,:) = GetThreeAxisByPlane(strG17, nxTempAxis, nLastThreeAxis);
    end
end

%%% Unified radius programming to center point programming
% Parameter strCurGId: G command (G2 or G3)
% Parameter nxLastPos: Start point of arc
% Parameter nxCurPos: End point of arc
% Parameter nRadius: Arc radius
% Return value nxArcCenter: Arc center
function nxArcCenter = TransRadiusToArcCenter(strCurGId, nxLastPos, nxCurPos, nRadius)
    % Calculate chord length
    nChordLength = norm(nxCurPos - nxLastPos);
    % Calculate the rotation angle
    nCosAngle = nChordLength * 0.5 / abs(nRadius);
    nSinAngle = sqrt(1 - nCosAngle * nCosAngle);
    % Calculate the direction of rotation
    if strcmp('G2', strCurGId) && nRadius > 0
        nSinAngle = -nSinAngle;
    end
    if strcmp('G3', strCurGId) && nRadius < 0
        nSinAngle = -nSinAngle;
    end
    % Calculate the center of arc
    nAlpha = abs(nRadius) / nChordLength;
    nxArcCenter = zeros(1,2);
    nxArcCenter(1) = nxLastPos(1) + nAlpha * ((nxCurPos(1) - nxLastPos(1)) * nCosAngle...
        - (nxCurPos(2) - nxLastPos(2)) * nSinAngle);
    nxArcCenter(2) = nxLastPos(2) + nAlpha * ((nxCurPos(1) - nxLastPos(1)) * nSinAngle...
        + (nxCurPos(2) - nxLastPos(2)) * nCosAngle);
end

%%% Convert 3D coordinate to 2D coordinate according to the coordinate plane
% Parameter strG17: Coordinate plane
% Parameter nxThreeAxis: Coordinate in three-dimension
% Return value nxTwoAxis: Coordinate in two-dimension
% Return value nThreeAxis: Coordinate in the third dimension
function [nxTwoAxis, nThreeAxis] = GetTwoAxisByPlane(strG17, nxThreeAxis)
    nxTwoAxis = zeros(1,2);
    switch strG17
        case 'G17'
            nxTwoAxis(1) = nxThreeAxis(1);
            nxTwoAxis(2) = nxThreeAxis(2);
            nThreeAxis = nxThreeAxis(3);
        case 'G18'
            nxTwoAxis(1) = nxThreeAxis(3);
            nxTwoAxis(2) = nxThreeAxis(1);
            nThreeAxis = nxThreeAxis(2);
        case 'G19'
            nxTwoAxis(1) = nxThreeAxis(2);
            nxTwoAxis(2) = nxThreeAxis(3);
            nThreeAxis = nxThreeAxis(1);
        otherwise
            error('Wrong G command: %s\n', strG17);
    end
end

%%% Convert 2D coordinate to 3D coordinate according to the coordinate plane
% Parameter strG17: Coordinate plane
% Parameter nxTwoAxis: Coordinate in two-dimension
% Parameter nThreeAxis: Coordinate in the third dimension
% Return value nxThreeAxis: Coordinate in three-dimension
function nxThreeAxis = GetThreeAxisByPlane(strG17, nxTwoAxis, nThreeAxis)
    nxThreeAxis = zeros(1,3);
    switch strG17
        case 'G17'
            nxThreeAxis(1) = nxTwoAxis(1);
            nxThreeAxis(2) = nxTwoAxis(2);
            nxThreeAxis(3) = nThreeAxis;
        case 'G18'
            nxThreeAxis(3) = nxTwoAxis(1);
            nxThreeAxis(1) = nxTwoAxis(2);
            nxThreeAxis(2) = nThreeAxis;
        case 'G19'
            nxThreeAxis(2) = nxTwoAxis(1);
            nxThreeAxis(3) = nxTwoAxis(2);
            nxThreeAxis(1) = nThreeAxis;
        otherwise
            error('Wrong G command: %s\n', strG17);
    end
end