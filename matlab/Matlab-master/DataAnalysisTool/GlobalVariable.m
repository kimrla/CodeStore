%% Global variables
function GlobalVariable(handles)
    % Parameters
    global g_period;            % Sampling period
    global g_dataFormat;        % Data format
    global g_validColumn;       % Valid data column
    global g_deleteRepeat;      % Delete redundant data at the beginning and end
    global g_SingleAxisMode;    % Single axis mode: In this mode, only single axis speed, acceleration, and jerk can be analyzed.
                                % Track and distance-speed cannot be drawn.
    global g_drawTrack;         % Draw track
    global g_drawTimeSpeed;     % Draw time-speed
    global g_drawTimeAcc;       % Draw time-acceleration
    global g_drawTimeJerk;      % Draw time-jerk
    global g_drawDistSpeed;     % Draw distance-speed
    % Initialize
    g_period = handles.metricdata.period;
    g_dataFormat = handles.metricdata.dataFormat;
    g_validColumn = textscan(handles.metricdata.validColumn, '%f');
    g_validColumn = g_validColumn{1};
    g_validColumn = sort(g_validColumn);
    g_deleteRepeat = handles.metricdata.deleteRepeat;
    g_SingleAxisMode = handles.metricdata.singleAxisMode;
    g_drawTrack = handles.metricdata.drawTrack;
    g_drawTimeSpeed = handles.metricdata.drawTimeSpeed;
    g_drawTimeAcc = handles.metricdata.drawTimeAcc;
    g_drawTimeJerk = handles.metricdata.drawTimeJerk;
    g_drawDistSpeed = handles.metricdata.drawDistSpeed;
end