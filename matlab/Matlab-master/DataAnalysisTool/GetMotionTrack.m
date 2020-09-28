%% Get motion track
function [motionTracks, motionFlags] = GetMotionTrack(allFiles)
    motionTracks = [];
    motionFlags = [];
    if isempty(allFiles)
        return;
    end
    % Analyse files
    if ~iscell(allFiles)
        [motionTrack, motionFlag] = GetTrack(allFiles);
        if ~isempty(motionTrack)
            motionTracks{1} = motionTrack;
            motionFlags{1} = motionFlag;
        end
    else
        motionTracks = cell(length(allFiles),1);
        motionFlags = cell(length(allFiles),1);
        index = 0;
        for i = 1:length(allFiles)
            [motionTrack, motionFlag] = GetTrack(allFiles{i});
            if ~isempty(motionTrack)
                index = index + 1;
                motionTracks{index} = motionTrack;
                motionFlags{index} = motionFlag;
            end
        end
        motionTracks = motionTracks(1:index);
        motionFlags = motionFlags(1:index);
    end
end

%%% Get motion track
function [motionTrack, motionFlag] = GetTrack(filePath)
    motionTrack = [];
    motionFlag = [];
    % Check if file exists
    fileExist = exist(filePath, 'file');
    if (0 == fileExist) || (7 == fileExist)
        fprintf('File (%s) does not exist!\n', filePath);
        return;
    end
    % file name
    fileInfo = dir(filePath);
    motionFlag = fileInfo.name;
    % track
    global g_dataFormat;
    fidRead = fopen(filePath, 'r');
    fAll = textscan(fidRead, g_dataFormat);
    fclose(fidRead);
    % Filter numeric columns
    validTrack = zeros(length(fAll{1}), length(fAll));
    index = 0;
    if size(validTrack,1) > 1
        for i = 1:length(fAll)
            if isnumeric(fAll{i}(1))
                index = index + 1;
                validTrack(:, index) = fAll{i};
            end
        end
    end
    validTrack = validTrack(:,1:index);
    % Filter for valid numeric column
    motionTrack = GetValidTrack(validTrack);
    % Delete redundant data at the beginning and end
    motionTrack = DeleteRepeatTrack(motionTrack);
end

%%% Filter for valid numeric column
function validTrack = GetValidTrack(origTrack)
    if size(origTrack,1) < 1
        validTrack = [];
        return;
    end
    global g_validColumn;
    validTrack = zeros(size(origTrack));
    index = 0;
    column = 1;
    for i = 1:size(origTrack,2)
        if column > length(g_validColumn)
            break;
        end
        if i == g_validColumn(column)
            index = index + 1;
            validTrack(:,index) = origTrack(:,i);
            column = column + 1;
        end
    end
    if index > 0
        validTrack = validTrack(:,1:index);
    else
        validTrack = [];
    end
end

%%% Delete redundant data at the beginning and end
function simplifyData = DeleteRepeatTrack(origData)
    global g_deleteRepeat;
    if ~g_deleteRepeat
        simplifyData = origData;
        return;
    end
    % Redundant data at the beginning
    startIndex = 1;
    while startIndex < size(origData,1)
        if ~isequal(origData(startIndex,:), origData(startIndex+1,:))
            break;
        end
        startIndex = startIndex + 1;
    end
    % Redundant data at the end
    endIndex = size(origData,1);
    while endIndex > 1
        if ~isequal(origData(endIndex,:), origData(endIndex-1,:))
            break;
        end
        endIndex = endIndex - 1;
    end
    % Cut data
    simplifyData = origData(startIndex:endIndex,:);
end