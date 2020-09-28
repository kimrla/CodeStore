%% Analyze motion track and draw figures based on options
function AnalysisMotionTrack(motionTracks, motionFlags)
    if isempty(motionTracks)
        fprintf('No motion track to analyze!\n');
        return;
    end
    % Drawing options
    global FileEntity;
    figId = FileEntity{end}.FigHandle;
    figColor = {'.-b', '.-g', '.-r', '.-m', '.-c', '.-y', '.-w'};
    % Drawing
    figLegend = '';
    for i = 1:min(length(motionTracks), length(figColor)-1)
        Id = DrawMotionTrack(motionTracks{i}, figId, figColor{i});
        figLegend = strcat(figLegend, '''', motionFlags{i}, '''', ',');
    end
    for i = figId+1:Id
        FileEntity{size(FileEntity,2)+1}.FigHandle = i;
        figure(i);
        eval(['legend(', figLegend(1:length(figLegend)-1), ')']);
    end
end

%%% Track drawing
function Id = DrawMotionTrack(motionTrack, figId, figColor)
    Id = figId;
    % track
    global g_drawTrack;
    if g_drawTrack
        Id = DrawTrack(motionTrack, Id, figColor);
    end
    % speed
    global g_drawTimeSpeed;
    global g_drawTimeAcc;
    global g_drawTimeJerk;
    global g_drawDistSpeed;
    if ~g_drawTimeSpeed && ~g_drawTimeAcc && ~g_drawTimeJerk && ~g_drawDistSpeed
        return;
    end
    motionSpeed = GetSpeed(motionTrack);
    if g_drawTimeSpeed
        Id = DrawMotionWithTitle(motionSpeed, Id, figColor, 'Time-Speed');
    end
    if g_drawDistSpeed
        Id = DrawDistSpeed(motionSpeed, Id, figColor);
    end
    % acceleration
    if ~g_drawTimeAcc && ~g_drawTimeJerk
        return;
    end
    motionAcc = GetAcc(motionSpeed);
    if g_drawTimeAcc
        Id = DrawMotionWithTitle(motionAcc, Id, figColor, 'Time-Acceleration');
    end
    % jerk
    if ~g_drawTimeJerk
        return;
    end
    motionJerk = GetAcc(motionAcc);
    Id = DrawMotionWithTitle(motionJerk, Id, figColor, 'Time-Jerk');
end

%%% Draw track
function Id = DrawTrack(motionTrack, figId, figColor)
    Id = figId;
    global g_SingleAxisMode;
    if g_SingleAxisMode
        fprintf('Track drawing is not supported in single axis mode!\n');
        return;
    end
    if size(motionTrack,2) < 1
        fprintf('No data for plotting!\n');
        return;
    end
    if size(motionTrack,2) > 3
        fprintf('Drawing of track exceeding 3 dimensions is not supported: the actual dimension is %d\n', size(motionTrack,2));
        return;
    end
    % Drawing
    Id = Id + 1;
    figure(Id);
    hold on;
    if size(motionTrack,2) == 1
        plot(motionTrack, zeros(size(motionTrack,1),1), figColor);
    elseif size(motionTrack,2) == 2
        plot(motionTrack(:,1), motionTrack(:,2), figColor);
    else
        plot3(motionTrack(:,1), motionTrack(:,2), motionTrack(:,3), figColor);
        view([1 1 1]);
    end
    title('Track');
    axis equal;
    grid on;
end

%%% Get speed
function motionSpeed = GetSpeed(motionTrack)
    global g_period;
    [row, column] = size(motionTrack);
    motionSpeed = zeros(row-1,column+2);
    for i = 1:row-1
        motionSpeed(i,1) = i * g_period;
        for j = 1:column
            motionSpeed(i,j+1) = (motionTrack(i+1,j) - motionTrack(i,j)) / g_period;
        end
        motionSpeed(i,column+2) = norm(motionSpeed(i,2:column+1));
    end
end

%%% Get acceleration
function motionAcc = GetAcc(motionSpeed)
    global g_period;
    [row, column] = size(motionSpeed);
    motionAcc = zeros(row-1,column);
    for i = 1:row-1
        motionAcc(i,1) = motionSpeed(i,1);
        for j = 2:column
            motionAcc(i,j) = (motionSpeed(i+1,j) - motionSpeed(i,j)) / g_period;
        end
    end
end

%%% Draw according to the title (speed, acceleration, jerk)
function Id = DrawMotionWithTitle(motionParam, figId, figColor, figTitle)
    Id = figId;
    global g_SingleAxisMode;
    % Determine figure layout based on dimensions: the first column is time, and the last column is the sum of all dimensions
    column = size(motionParam,2);
    if column < 3
        fprintf('Wrong parameter: insufficient data dimensions!\n');
        return;
    end
    preTitle = {'X axis ', 'Y axis ', 'Z axis ', 'Combined '};
    if column == 3
        % One dimensional data: only draw single axis data
        Id = Id + 1;
        figure(Id);
        hold on;
        plot(motionParam(:,1), motionParam(:,2), figColor);
        title([preTitle{1}, figTitle]);
        axis auto; grid on;
    elseif column == 4 && ~g_SingleAxisMode
        % Two dimensional data: draw single axis data and combined data
        Id = Id + 1;
        figure(Id);
        handleSub = cell(3,1);
        for i = 1:3
            handleSub{i} = subplot(3,1,i);
            hold on;
            plot(motionParam(:,1), motionParam(:,i+1), figColor);
            if i ~= 3
                j = i;
            else
                j = i + 1;
            end
            title([preTitle{j}, figTitle]);
            axis auto; grid on;
        end
        linkaxes([handleSub{1}, handleSub{2}, handleSub{3}], 'x');
    elseif column == 5 && ~g_SingleAxisMode
        % Three dimensional data: draw single axis data and combined data
        for i = 1:2:3
            Id = Id + 1;
            figure(Id);
            handleSub = cell(2,1);
            for j = 1:2
                handleSub{j} = subplot(2,1,j);
                hold on;
                plot(motionParam(:,1), motionParam(:,i+j), figColor);
                title([preTitle{i+j-1}, figTitle]);
                axis auto; grid on;
            end
            linkaxes([handleSub{1}, handleSub{2}], 'x');
        end
    else
        % Draw only single axis data: data exceeding three dimensions does not draw data by default
        for i = 1:2:column-3
            Id = Id + 1;
            figure(Id);
            handleSub = cell(2,1);
            for j = 1:2
                handleSub{j} = subplot(2,1,j);
                hold on;
                plot(motionParam(:,1), motionParam(:,i+j), figColor);
                title(['No. ', num2str(i+j-1), ' ', figTitle]);
                axis auto; grid on;
            end
            linkaxes([handleSub{1}, handleSub{2}], 'x');
        end
        if mod(column,2)
            Id = Id + 1;
            figure(Id);
            hold on;
            plot(motionParam(:,1), motionParam(:,column-1), figColor);
            title(['No. ', num2str(column-1), ' ', figTitle]);
            axis auto; grid on;
        end
    end
end

%%% Draw distance-speed
function Id = DrawDistSpeed(motionSpeed, figId, figColor)
    Id = figId;
    global g_SingleAxisMode;
    if g_SingleAxisMode
        fprintf('Distance-speed drawing is not supported in single axis mode!\n');
        return;
    end
    [row, column] = size(motionSpeed);
    if column < 3 || column > 5
        error('%d dimension data does not support drawing distance-speed!\n', column-2);
    end
    global g_period;
    distSpeed = zeros(row,2);
    distSpeed(1,1) = motionSpeed(1,column) * g_period;
    distSpeed(1,2) = motionSpeed(1,column);
    for i = 2:row
        distSpeed(i,1) = distSpeed(i-1,1) + motionSpeed(i,column) * g_period;
        distSpeed(i,2) = motionSpeed(i,column);
    end
    % Drawing
    Id = figId + 1;
    figure(Id);
    hold on;
    plot(distSpeed(:,1), distSpeed(:,2), figColor);
    title('Distance-Speed');
    grid on;
end