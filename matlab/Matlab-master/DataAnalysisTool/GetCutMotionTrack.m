%% Get cut motion track
function [cutMotionTracks, cutMotionFlags] = GetCutMotionTrack(filePath, cutIndexs)
    cutMotionTracks = [];
    cutMotionFlags = [];
    [motionTracks, motionFlags] = GetMotionTrack(filePath);
    % Parameter judgment
    if length(motionTracks) ~= 1 || length(motionTracks) ~= length(motionFlags)
        fprintf('The motion track does not meet the cutting conditions: the number of motion track is %d\n', length(motionTracks));
        return;
    end
    motionTrack = motionTracks{1};
    motionFlag = motionFlags{1};
    % Get cut track
    global g_period;
    cutMotionTracks = cell(length(cutIndexs),1);
    cutMotionFlags = cell(length(cutIndexs),1);
    index = 0;
    for i = 1:length(cutIndexs)
        startIndex = max(cutIndexs{i}.startIndex,1);
        endIndex = min(cutIndexs{i}.endIndex, size(motionTrack,1));
        if startIndex >= endIndex
            continue;
        end
        index = index + 1;
        cutMotionTracks{index} = motionTrack(startIndex:endIndex,:);
        cutMotionFlags{index} = [motionFlag, ': ', num2str(cutIndexs{i}.startIndex * g_period), ...
            '-', num2str(cutIndexs{i}.endIndex * g_period)];
    end
    if index > 0
        cutMotionTracks = cutMotionTracks(1:index);
        cutMotionFlags = cutMotionFlags(1:index);
    else
        cutMotionTracks = [];
        cutMotionFlags = [];
    end
end