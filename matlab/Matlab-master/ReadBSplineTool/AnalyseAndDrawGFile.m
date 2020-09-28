%% Analyze the motion instructions in the tool path and draw pictures
% Parameter filePath: Toolpath file path
% Parameter nFigId: Picture number
% Parameter nColorId: Color number
% Return value figHandle: Figure handle
function figHandle = AnalyseAndDrawGFile(filePath, nFigId, nColorId)
    global g_nHandleCount;
    global g_figHandle;
    g_nHandleCount = 4;
    g_figHandle = zeros(g_nHandleCount,1);
    % Open a file
    fidRead = fopen(filePath, 'r');
    % Programming mode: absolute value programming G90, relative value programming G91
    strG90 = '';
    % Coordinate plane: XY plane G17, ZX plane G18, YZ plane G19
    strG17 = '';
    % Motion G command identifier: G0, G1, G2, G3, G6.2
    strLastGId = '';
    strCurGId = '';
    % Coordinate values: x-axis coordinate, y-axis coordinate, z-axis coordinate
    bSetStartPos = false;
    nxLastPos = ones(1,3) * Inf;
    nxCurPos = ones(1,3) * Inf;
    % Parameters: circle center I, circle center J, circle center / nurbs node knot K, arc radius / nurbs node weight R, nurbs node degree P
    nxCurParam = ones(1,5) * Inf;
    % Parameters at the end of the program early
    bForceEnd = false;
    % B-spline node: degree, knot vector, poles and weights
    nurbs.nDegree = 0;
    nurbs.vecKnots = [];
    nurbs.vecPoles = [];
    nurbs.vecWeights = [];
    % Cached data points
    global g_nMaxCachePoints;
    g_nMaxCachePoints = 10000;
    InitCellPointCache(g_nMaxCachePoints, nColorId);
    % Curve scatter precision
    global g_nScatterPrecision;
    % Read files in batches until the end of the file is reached
    while ~feof(fidRead)
        % Read multiple lines of file at once
        fAll = textscan(fidRead, '%s', 1000, 'Delimiter', '\n');
        fCell = fAll{1};
        nMaxLines = size(fCell,1);
        for i = 1:nMaxLines
            % Get and process each line of string
            strOrig = fCell{i};
            % Judge whether if it is a blank line
            if isempty(strOrig)
                continue;
            end
            % Judge whether if G command is included
            if isempty(strLastGId) && ~contains(strOrig, 'G')
                % If the tool path contains coordinates, the G command defaults to G0
                if contains(strOrig, 'X') || contains(strOrig, 'Y') || contains(strOrig, 'Z')
                    strCurGId = 'G0';
                else
                    continue;
                end
            end
            % Judge whether if the program is over
            bForceEnd = JudgeEndMachine(strOrig);
            if bForceEnd
                break;
            end
            % Read G command and coordinate value
            [strG90, strG17, strCurGId, nxCurPos, nxCurParam] =...
                AnalyseGCode(strOrig, strG90, strG17, strCurGId, nxCurPos, nxCurParam);
            % Judge the type of G command
            if isempty(strLastGId) && isempty(strCurGId)
                continue;
            elseif isempty(strCurGId)
                strCurGId = strLastGId;
            end
            % Convert invalid coordinates to 0
            if isinf(nxCurPos(1)) && isinf(nxCurPos(2)) && isinf(nxCurPos(3)) && isinf(nxCurParam(3))
                if ~isempty(strCurGId)
                    strLastGId = strCurGId;
                end
                continue;
            else
                if bSetStartPos
                    nxCurPos = GetAbsoluteAxis(strG90, nxLastPos, nxCurPos);
                else
                    nxLastPos = GetAbsoluteAxis(strG90, zeros(1,3), nxCurPos);
                    nxCurPos = nxLastPos;
                    bSetStartPos = true;
                end
            end
            % Check if the previous command is a B-spline node
            if ~strcmp(strLastGId, strCurGId) && strcmp(strLastGId, 'G6.2')
                % Scatter nurbs node
                nurbs.vecPoles = nurbs.vecPoles(1:size(nurbs.vecPoles,1)-nurbs.nDegree-1,:);
                vecPoints = ScatterNurbs(nurbs, g_nScatterPrecision);
                StoreCellPointCache(vecPoints, 3, nFigId);
                StoreCellPointCache(nurbs.vecPoles, 4, nFigId);
                nurbs.nDegree = 0;
                nurbs.vecKnots = [];
                nurbs.vecPoles = [];
                nurbs.vecWeights = [];
            end
            % Analyze coordinate values
            if strcmp(strCurGId, 'G1')
                StoreCellPointCache([nxLastPos; nxCurPos], 1, nFigId);
            elseif strcmp(strCurGId, 'G2') || strcmp(strCurGId, 'G3')
                vecPoints = ScatterArc(strG17, strCurGId, nxLastPos, nxCurPos, nxCurParam, g_nScatterPrecision);
                StoreCellPointCache(vecPoints, 2, nFigId);
            elseif strcmp(strCurGId, 'G6.2')
                % The degree of nurbs node
                if ~isinf(nxCurParam(5))
                    % Check if there are two connected B-spline nodes
                    if strcmp(strLastGId, 'G6.2')
                        % Scatter nurbs node
                        nurbs.vecPoles = nurbs.vecPoles(1:size(nurbs.vecPoles,1)-nurbs.nDegree-1,:);
                        vecPoints = ScatterNurbs(nurbs, g_nScatterPrecision);
                        StoreCellPointCache(vecPoints, 3, nFigId);
                        StoreCellPointCache(nurbs.vecPoles, 4, nFigId);
                        nurbs.vecKnots = [];
                        nurbs.vecPoles = [];
                        nurbs.vecWeights = [];
                    end
                    nurbs.nDegree = nxCurParam(5);
                end
                % The knot vector of nurbs node
                if ~isinf(nxCurParam(3))
                    nurbs.vecKnots = [nurbs.vecKnots; nxCurParam(3)];
                end
                % The poles of nurbs node
                if 0 == size(nurbs.vecPoles,1)
                    nurbs.vecPoles = nxCurPos;
                    % The weights of poles
                    if ~isinf(nxCurParam(4))
                        nurbs.vecWeights = nxCurParam(4);
                    end
                else
                    nurbs.vecPoles = [nurbs.vecPoles; nxCurPos];
                    if ~isinf(nxCurParam(4))
                        nurbs.vecWeights = [nurbs.vecWeights; nxCurParam(4)];
                    end
                end
            else
                if ~strcmp(strCurGId, 'G0')
                    error('Wrong G command: %s\n', strCurGId);
                end
            end
            % Record last parameter value
            strLastGId = strCurGId;
            nxLastPos = nxCurPos;
            % Parameters reset
            nxCurPos = ones(1,3) * Inf;
            nxCurParam = ones(1,5) * Inf;
        end
        if bForceEnd
            break;
        end
    end
    % The last paragraph is a B-spline node
    if strcmp(strLastGId, 'G6.2')
        % Scatter nurbs node
        nurbs.vecPoles = nurbs.vecPoles(1:size(nurbs.vecPoles,1)-nurbs.nDegree-1,:);
        vecPoints = ScatterNurbs(nurbs, g_nScatterPrecision);
        StoreCellPointCache(vecPoints, 3, nFigId);
        StoreCellPointCache(nurbs.vecPoles, 4, nFigId);
    end
    % Flush cache
    for i = 1:g_nHandleCount
        FlushCellPointCache(i, nFigId);
    end
    % Close file
    fclose(fidRead);
    figHandle = g_figHandle;
end

%%% Initialize cache data points
function InitCellPointCache(nMaxSize, nColorId)
    % Check parameters
    if nColorId ~= 1 && nColorId ~= 2
        error('Wrong parameter: color number is %f\n', nColorId);
    end
    global g_vecColor;
    global g_nHandleCount;
    global g_cellPointCache;
    g_cellPointCache = cell(g_nHandleCount,1);
    for i = 1:g_nHandleCount
        % Data point
        g_cellPointCache{i}.points = zeros(nMaxSize,3);
        % Start index
        g_cellPointCache{i}.index = 0;
        % Figure color
        g_cellPointCache{i}.color = g_vecColor{nColorId}{i};
    end
    global g_nxMaxRange;
    global g_nxMinRange;
    global g_bFirstPoint;
    g_nxMaxRange = zeros(1,3);
    g_nxMinRange = zeros(1,3);
    g_bFirstPoint = 0;
end

%%% Store cache data points
function StoreCellPointCache(vecPoints, nHandleIndex, nFigId)
    % Check parameters
    if size(vecPoints,1) < 1
        error('Wrong parameter: number of data points is %f\n', size(vecPoints,1));
    end
    global g_nHandleCount;
    if nHandleIndex < 1 || nHandleIndex > g_nHandleCount
        error('Wrong parameter: handle index is %f\n', nHandleIndex);
    end
    % If the colors are consistent, modify the handle index value
    global g_cellPointCache;
    for i = 1:g_nHandleCount-1
        if strcmp(g_cellPointCache{i}.color, g_cellPointCache{nHandleIndex}.color)
            nHandleIndex = i;
            break;
        end
    end
    % Store data points
    nPointIndex = g_cellPointCache{nHandleIndex}.index;
    if nPointIndex > 0 && norm(g_cellPointCache{nHandleIndex}.points(nPointIndex,:) - vecPoints(1,:)) > 0.001
        % The current data point and the cached data point are not connected end to end
        vecNewPoints = g_cellPointCache{nHandleIndex}.points(1:nPointIndex,:);
        PlotFigureByPlane(vecNewPoints, nHandleIndex, nFigId, g_cellPointCache{nHandleIndex}.color, 0);
        g_cellPointCache{nHandleIndex}.index = 0;
    end
    % Check if the cache can fit
    global g_nMaxCachePoints;
    nPointIndex = g_cellPointCache{nHandleIndex}.index;
    if nPointIndex + size(vecPoints,1) >= g_nMaxCachePoints
        vecNewPoints = [g_cellPointCache{nHandleIndex}.points(1:nPointIndex,:); vecPoints];
        PlotFigureByPlane(vecNewPoints, nHandleIndex, nFigId, g_cellPointCache{nHandleIndex}.color, 0);
        g_cellPointCache{nHandleIndex}.index = 0;
        return;
    end
    % Cache data points
    if nPointIndex == 0
        g_cellPointCache{nHandleIndex}.points = vecPoints;
        g_cellPointCache{nHandleIndex}.index = size(vecPoints,1);
    else
        g_cellPointCache{nHandleIndex}.points(nPointIndex+1:nPointIndex+size(vecPoints,1)-1,:) = vecPoints(2:end,:);
        g_cellPointCache{nHandleIndex}.index = g_cellPointCache{nHandleIndex}.index + size(vecPoints,1) - 1;
    end
end

%%% Flush cache data points
function FlushCellPointCache(nHandleIndex, nFigId)
    % Check parameters
    global g_nHandleCount;
    if nHandleIndex < 1 || nHandleIndex > g_nHandleCount
        error('Wrong parameter: handle index is %f\n', nHandleIndex);
    end
    global g_cellPointCache;
    nPointIndex = g_cellPointCache{nHandleIndex}.index;
    if nPointIndex == 0
        return;
    end
    vecNewPoints = g_cellPointCache{nHandleIndex}.points(1:nPointIndex,:);
    PlotFigureByPlane(vecNewPoints, nHandleIndex, nFigId, g_cellPointCache{nHandleIndex}.color, 1);
    g_cellPointCache{nHandleIndex}.index = 0;
end

%%% Judge whether the program ends processing
% Parameter strOrig: String to parse
% Return value bForceEnd: Whether to stop processing, true to stop, false to not stop
function bForceEnd = JudgeEndMachine(strOrig)
    bForceEnd = false;
    % If M2 or M30 is included, the program ends
    nIndex = strfind(strOrig, 'M');
    if isempty(nIndex)
        nIndex = strfind(strOrig, 'm');
        if isempty(nIndex)
            return;
        end
    end
    % The first character is the identifier by default
    nResult = textscan(strOrig(nIndex:end), '%c%f', 1);
    if ~isempty(nResult{2})
        nValue = nResult{2};
        if (2 == nValue) || (30 == nValue)
            bForceEnd = true;
        end
    else
        error('Index value does not exist: string is %s, index is %d\n', strOrig, nIndex);
    end
end

%%% Get absolute coordinate
% Parameter strG90: Programming mode, G90 is absolute coordinate programming (default value), G91 is relative coordinate programming
% Parameter nxLastPos: Absolute coordinates of the previous interpolation point
% Parameter nxCurPos: Absolute or relative coordinates of the current interpolation point
% Return value nxCurPosOut: Get the absolute coordinates of the current interpolation point according to the programming mode
function nxCurPosOut = GetAbsoluteAxis(strG90, nxLastPos, nxCurPos)
    % Check parameters
    if ~isempty(strG90) && strcmp(strG90, 'G90') && strcmp(strG90, 'G91')
        error('Wrong programming: %s\n', strG90);
    end
    if isinf(nxLastPos(1)) && isinf(nxLastPos(2)) && isinf(nxLastPos(3))
        nxLastPos = zeros(1,3);
    end
    if isinf(nxLastPos(1)) || isinf(nxLastPos(2)) || isinf(nxLastPos(3))
        error('Wrong coordinate: X = %f, Y = %f, Z = %f\n', nxLastPos(1), nxLastPos(2), nxLastPos(3));
    end
    % Absolute coordinate
    nxCurPosOut = nxLastPos;
    for i = 1:3
        if ~isinf(nxCurPos(i)) && strcmp(strG90, 'G91')
            nxCurPosOut(i) = nxLastPos(i) + nxCurPos(i);
        elseif ~isinf(nxCurPos(i)) && ~strcmp(strG90, 'G91')
            nxCurPosOut(i) = nxCurPos(i);
        end
    end
end

%%% Analyze G commands and coordinate values
% Parameter strOrig: String to parse
% Parameter strG90: Programming style
% Parameter strG17: Coordinate plane
% Parameter strCurGId: G command identifier
% Parameter nxCurPos: Coordinate value
% Parameter nxCurParam: Arc / B-spline parameters
% Return value: The outout parameters are similar to the input parameters
function [strG90Out, strG17Out, strCurGIdOut, nxCurPosOut, nxCurParamOut]...
    = AnalyseGCode(strOrig, strG90, strG17, strCurGId, nxCurPos, nxCurParam)
    % Defaults
    strG90Out = strG90;
    strG17Out = strG17;
    strCurGIdOut = strCurGId;
    nxCurPosOut = nxCurPos;
    nxCurParamOut = nxCurParam;
    % Parsing values in the string
    vecResult = textscan(strOrig, '%c%f');
    vecFlags = vecResult{1};
    vecValues = vecResult{2};
    for i = 1:length(vecFlags)
        if strcmp(vecFlags(i), 'G') || strcmp(vecFlags(i), 'g')
            % Classification processing G instruction identifier
            switch vecValues(i)
                case 1
                    strCurGIdOut = 'G1';
                case 2
                    strCurGIdOut = 'G2';
                case 3
                    strCurGIdOut = 'G3';
                case 6.2
                    strCurGIdOut = 'G6.2';
                case 0
                    strCurGIdOut = 'G0';
                case 17
                    strG17Out = 'G17';
                case 18
                    strG17Out = 'G18';
                case 19
                    strG17Out = 'G19';
                case 90
                    strG90Out = 'G90';
                case 91
                    strG90Out = 'G91';
            end
        elseif strcmp(vecFlags(i), 'X') || strcmp(vecFlags(i), 'x')
            nxCurPosOut(1) = vecValues(i);
        elseif strcmp(vecFlags(i), 'Y') || strcmp(vecFlags(i), 'y')
            nxCurPosOut(2) = vecValues(i);
        elseif strcmp(vecFlags(i), 'Z') || strcmp(vecFlags(i), 'z')
            nxCurPosOut(3) = vecValues(i);
        elseif strcmp(vecFlags(i), 'I') || strcmp(vecFlags(i), 'i')
            nxCurParamOut(1) = vecValues(i);
        elseif strcmp(vecFlags(i), 'J') || strcmp(vecFlags(i), 'j')
            nxCurParamOut(2) = vecValues(i);
        elseif strcmp(vecFlags(i), 'K') || strcmp(vecFlags(i), 'k')
            nxCurParamOut(3) = vecValues(i);
        elseif strcmp(vecFlags(i), 'R') || strcmp(vecFlags(i), 'r')
            nxCurParamOut(4) = vecValues(i);
        elseif strcmp(vecFlags(i), 'P') || strcmp(vecFlags(i), 'p')
            nxCurParamOut(5) = vecValues(i);
        end
    end
end

%%% Drawing based on coordinate view
function PlotFigureByPlane(vecFitPoint, nHandleIndex, nFigId, strColor, bFlush)
    % Check parameters
    global g_nHandleCount;
    if nHandleIndex < 1 || nHandleIndex > g_nHandleCount
        error('Wrong parameter: handle index is %f\n', nHandleIndex);
    end
    global g_bDrawControlPoint;
    if nHandleIndex == g_nHandleCount
        if ~g_bDrawControlPoint
            return;
        end
    end
    % Record maximum and minimum range values
    global g_nxMaxRange;
    global g_nxMinRange;
    global g_bFirstPoint;
    nxMaxRange = [max(vecFitPoint(:,1)), max(vecFitPoint(:,2)), max(vecFitPoint(:,3))];
    nxMinRange = [min(vecFitPoint(:,1)), min(vecFitPoint(:,2)), min(vecFitPoint(:,3))];
    if ~g_bFirstPoint
        g_bFirstPoint = 1;
        g_nxMaxRange = nxMaxRange;
        g_nxMinRange = nxMinRange;
    else
        for i = 1:3
            g_nxMaxRange(i) = max([g_nxMaxRange(i), nxMaxRange(i)]);
            g_nxMinRange(i) = min([g_nxMinRange(i), nxMinRange(i)]);
        end
    end
    if bFlush
        nSpace = max([abs(g_nxMaxRange(1) - g_nxMinRange(1)), abs(g_nxMaxRange(2) - g_nxMinRange(2)),...
            abs(g_nxMaxRange(3) - g_nxMinRange(3))]) * 0.05;
    end
    % Drawing
    global g_figHandle;
    global g_strPlane;
    figure(nFigId);
    hold on;
    switch g_strPlane
        case '3D View'
            if bFlush
                plot3(g_nxMaxRange(:,1)+nSpace, g_nxMaxRange(:,2)+nSpace, g_nxMaxRange(:,3)+nSpace, 'w');
                plot3(g_nxMinRange(:,1)-nSpace, g_nxMinRange(:,2)-nSpace, g_nxMinRange(:,3)-nSpace, 'w');
            end
            g_figHandle(nHandleIndex) = plot3(vecFitPoint(:,1), vecFitPoint(:,2), vecFitPoint(:,3), strColor);
            view([1 1 1]);
        case 'X-Y Plane'
            if bFlush
                plot(g_nxMaxRange(:,1)+nSpace, g_nxMaxRange(:,2)+nSpace, 'w');
                plot(g_nxMinRange(:,1)-nSpace, g_nxMinRange(:,2)-nSpace, 'w');
            end
            g_figHandle(nHandleIndex) = plot(vecFitPoint(:,1), vecFitPoint(:,2), strColor);
        case 'Y-Z Plane'
            if bFlush
                plot(g_nxMaxRange(:,2)+nSpace, g_nxMaxRange(:,3)+nSpace, 'w');
                plot(g_nxMinRange(:,2)-nSpace, g_nxMinRange(:,3)-nSpace, 'w');
            end
            g_figHandle(nHandleIndex) = plot(vecFitPoint(:,2), vecFitPoint(:,3), strColor);
        case 'Z-X Plane'
            if bFlush
                plot(g_nxMaxRange(:,3)+nSpace, g_nxMaxRange(:,1)+nSpace, 'w');
                plot(g_nxMinRange(:,3)-nSpace, g_nxMinRange(:,1)-nSpace, 'w');
            end
            g_figHandle(nHandleIndex) = plot(vecFitPoint(:,3), vecFitPoint(:,1), strColor);
        otherwise
            error('Wrong view: %s\n', g_strPlane);
    end
    grid on;
    axis equal;
end