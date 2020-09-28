%% test for nurbs node
% The nurbs node information is as follows:
%   nurbs.nDegree ------ The degree of nurbs node
%   nurbs.vecKnots ----- The knot vector of nurbs node
%   nurbs.vecPoles ----- The poles of nurbs node
%   nurbs.vecWeights --- The weights of nurbs node
function Test_Nurbs()
    clc; close all;
    % Read nurbs node from a file
    nurbs = ReadNurbsFromFile('TestData\circle.nc');
    if isempty(nurbs)
        return;
    end
    % Test scatter for nurbs node
    nCount = 200;
    nDeflection = 1e-3;
    Test_NurbsScatter(nurbs, nCount, nDeflection);
    % Test derivative for nurbs node
    Test_NurbsDeriv(nurbs, nCount);
end

%%% Read nurbs node from a file
function nurbs = ReadNurbsFromFile(strFilePath)
    nurbs = [];
    % Check if file exists
    bFileExist = exist(strFilePath, 'file');
    if 0 == bFileExist || 7 == bFileExist
        fprintf('File (%s) not exists!\n', strFilePath);
        return;
    end
    % Nurbs node information
    nurbs.nDegree = 0;
    nurbs.vecKnots = [];
    nurbs.vecPoles = [];
    nurbs.vecWeights = [];
    % Read file
    fidRead = fopen(strFilePath, 'r');
    while ~feof(fidRead)
        % Read multiple rows at once
        fAll = textscan(fidRead, '%s', 1000, 'Delimiter', '\n');
        fCell = fAll{1};
        nMaxLines = size(fCell,1);
        for i = 1:nMaxLines
            % Get and process each line of string
            strOrig = fCell{i};
            if isempty(strOrig)
                continue;
            end
            nurbs = GetNurbsFlags(nurbs, strOrig);
        end
    end
    fclose(fidRead);
    % Check if nurbs node is legal
    CheckNurbs(nurbs);
end

%%% Get the value of the nurbs node identifier
function nurbsNew = GetNurbsFlags(nurbs, strOrig)
    nurbsNew = nurbs;
    % Coordinate value
    nxPole = [0, 0, 0];
    bHasPole = false;
    % Parse identifiers and corresponding values
    vecResult = textscan(strOrig, '%c%f');
    vecFlags = vecResult{1};
    vecValues = vecResult{2};
    for i = 1:length(vecFlags)
        switch vecFlags(i)
            case 'G'
                if vecValues(i) ~= 6.2
                    error('Wrong parameter: G code identified is %f\n', vecValues(i));
                end
            case 'P'
                if nurbsNew.nDegree ~= 0
                    error('Only supports the analysis of single-segment nurbs node!\n');
                end
                nurbsNew.nDegree = vecValues(i);
            case 'K'
                nurbsNew.vecKnots = [nurbsNew.vecKnots; vecValues(i)];
            case 'X'
                nxPole(1) = vecValues(i);
                bHasPole = true;
            case 'Y'
                nxPole(2) = vecValues(i);
                bHasPole = true;
            case 'Z'
                nxPole(3) = vecValues(i);
                bHasPole = true;
            case 'R'
                nurbsNew.vecWeights = [nurbsNew.vecWeights; vecValues(i)];
            otherwise
                error('Wrong identifier : %c\n', vecFlags(i));
        end
    end
    if bHasPole
        nurbsNew.vecPoles = [nurbsNew.vecPoles; nxPole];
    end
end

%%% Test scatter for nurbs node
function Test_NurbsScatter(nurbs, nCount, nDeflection)
    % Calculate points on the curve
    vecPointsByFormula = zeros(nCount+1,3);
    nStartKnot = nurbs.vecKnots(nurbs.nDegree+1);
    nEndKnot = nurbs.vecKnots(length(nurbs.vecKnots)-nurbs.nDegree);
    for i = 0:nCount
        nKnot = nStartKnot + (nEndKnot - nStartKnot) * i / nCount;
        if size(nurbs.vecPoles,2) == 2
            [vecPointsByFormula(i+1,1:2), ~, ~] = GetNurbsDeriv(nurbs, nKnot);
        elseif size(nurbs.vecPoles,2) == 3
            [vecPointsByFormula(i+1,:), ~, ~] = GetNurbsDeriv(nurbs, nKnot);
        end
    end
    vecPointsByDeflection = ScatterNurbs(nurbs, nDeflection);
    % Drawing comparison
    figure();
    hold on;
    if isequal(vecPointsByFormula(:,3), zeros(size(vecPointsByFormula,1),1))
        plot(vecPointsByFormula(:,1), vecPointsByFormula(:,2), '.-b');
        plot(vecPointsByDeflection(:,1), vecPointsByDeflection(:,2), '.-g');
        plot(nurbs.vecPoles(:,1), nurbs.vecPoles(:,2), '*-r');
    else
        plot3(vecPointsByFormula(:,1), vecPointsByFormula(:,2), vecPointsByFormula(:,3), '.-b');
        plot3(vecPointsByDeflection(:,1), vecPointsByDeflection(:,2), vecPointsByDeflection(:,3), '.-g');
        plot3(nurbs.vecPoles(:,1), nurbs.vecPoles(:,2), nurbs.vecPoles(:,3), '*-r');
        view([1 1 1]);
    end
    legend('Scatting by uniformly sampled node vector', 'Scatting by bow height error', 'Poles of nurbs node');
    title('Track figure');
    grid on;
    axis equal;
end

%%% Test derivative for nurbs node
function Test_NurbsDeriv(nurbs, nCount)
    % Calculate derivative value
    vecNewDeriv = zeros(nCount+1,9);
    nStartKnot = nurbs.vecKnots(nurbs.nDegree+1);
    nEndKnot = nurbs.vecKnots(length(nurbs.vecKnots)-nurbs.nDegree);
    for i = 0:nCount
        nKnot = nStartKnot + (nEndKnot - nStartKnot) * i / nCount;
        if size(nurbs.vecPoles,2) == 2
            [vecNewDeriv(i+1,1:2), vecNewDeriv(i+1,4:5), vecNewDeriv(i+1,7:8)] = GetNurbsDeriv(nurbs, nKnot);
        elseif size(nurbs.vecPoles,2) == 3
            [vecNewDeriv(i+1,1:3), vecNewDeriv(i+1,4:6), vecNewDeriv(i+1,7:9)] = GetNurbsDeriv(nurbs, nKnot);
        end
    end
    nStepKnot = (nEndKnot - nStartKnot) / nCount;
    vecOrigDeriv = vecNewDeriv;
    vecOrigDeriv(2:end,4:6) = diff(vecOrigDeriv(:,1:3)) / nStepKnot;
    vecOrigDeriv(3:end,7:9) = diff(vecOrigDeriv(2:end,4:6)) / nStepKnot;
    % Drawing comparison
    vecTitles = {'First derivative figure', 'Second derivative figure'};
    for i = 2:min(nurbs.nDegree+1,3)
        figure();
        hold on;
        if isequal(vecNewDeriv(:,i*3), zeros(size(vecNewDeriv,1),1))
            plot(vecOrigDeriv(i:end,i*3-2), vecOrigDeriv(i:end,i*3-1), '.-b');
            plot(vecNewDeriv(:,i*3-2), vecNewDeriv(:,i*3-1), '.-g');
        else
            plot3(vecOrigDeriv(i:end,i*3-2), vecOrigDeriv(i:end,i*3-1), vecOrigDeriv(i:end,i*3), '.-b');
            plot3(vecNewDeriv(:,i*3-2), vecNewDeriv(:,i*3-1), vecNewDeriv(:,i*3), '.-g');
            view([1 1 1]);
        end
        legend('Differential method', 'Calculation method');
        title(vecTitles{i-1});
        grid on;
        axis equal;
    end
    % Calculating curvature
    vecNewCurvature = zeros(size(vecNewDeriv,1),1);
    vecOldCurvature = zeros(size(vecOrigDeriv,1),1);
    for i = 1:length(vecNewCurvature)
        nxNumerator = zeros(1,3);
        nxNumerator(1) = vecNewDeriv(i,5) * vecNewDeriv(i,9) - vecNewDeriv(i,8) * vecNewDeriv(i,6);
        nxNumerator(2) = vecNewDeriv(i,6) * vecNewDeriv(i,7) - vecNewDeriv(i,9) * vecNewDeriv(i,4);
        nxNumerator(3) = vecNewDeriv(i,4) * vecNewDeriv(i,8) - vecNewDeriv(i,7) * vecNewDeriv(i,5);
        vecNewCurvature(i) = norm(nxNumerator) / (norm(vecNewDeriv(i,4:6)))^3;
    end
    for i = 1:length(vecOldCurvature)
        nxNumerator = zeros(1,3);
        nxNumerator(1) = vecOrigDeriv(i,5) * vecOrigDeriv(i,9) - vecOrigDeriv(i,8) * vecOrigDeriv(i,6);
        nxNumerator(2) = vecOrigDeriv(i,6) * vecOrigDeriv(i,7) - vecOrigDeriv(i,9) * vecOrigDeriv(i,4);
        nxNumerator(3) = vecOrigDeriv(i,4) * vecOrigDeriv(i,8) - vecOrigDeriv(i,7) * vecOrigDeriv(i,5);
        vecOldCurvature(i) = norm(nxNumerator) / (norm(vecOrigDeriv(i,4:6)))^3;
    end
    % Drawing
    figure();
    hold on;
    plot(1:length(vecOldCurvature), vecOldCurvature, '.-b');
    plot(1:length(vecNewCurvature), vecNewCurvature, '.-g');
    legend('Differential method', 'Calculation method');
    title('Curvature figure');
    grid on;
end