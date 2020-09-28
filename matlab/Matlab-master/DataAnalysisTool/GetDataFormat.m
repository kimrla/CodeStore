%% Get data format
% Assumes the same number of columns in each row of data, but there is no limit to the number of columns
function strDataFormat = GetDataFormat(strFilePath)
    strDataFormat = [];
    % Read data file and parse data format
    fidRead = fopen(strFilePath, 'r');
    fAll = textscan(fidRead, '%s', 10, 'Delimiter', '\n');
    fclose(fidRead);
    fCell = fAll{1};
    nMaxLines = size(fCell,1);
    for i = 1:nMaxLines
        strOrig = fCell{i};
        if isempty(strOrig)
            continue;
        end
        strDataFormat = GetDataFormatFromString(strOrig);
        if ~isempty(strDataFormat)
            break;
        end
    end
end

%%% Get data format: parse the number of columns from a string
function strDataFormat = GetDataFormatFromString(strOrigData)
    % Assembly data format. The default separators are spaces, tabs, commas, and semicolons.
    strOrigData = strtrim(strOrigData);
    strDataFormat = '';
    for i = 1:length(strOrigData)
        if strcmp(strOrigData(i), char(32)) || strcmp(strOrigData(i), char(9))
            strDataFormat = strcat(strDataFormat, 32);
        elseif strcmp(strOrigData(i), ',')
            strDataFormat = strcat(strDataFormat, ',');
        elseif strcmp(strOrigData(i), ';')
            strDataFormat = strcat(strDataFormat, ';');
        else
            if isempty(strDataFormat)
                strDataFormat = strcat(strDataFormat, '%f');
            elseif ~strcmp(strDataFormat(end-1:end), '%f')
                if strcmp(strDataFormat(end), ' ')
                    strDataFormat = strcat(strDataFormat, 32, '%f');
                else
                    strDataFormat = strcat(strDataFormat, '%f');
                end
            end
        end
    end
end