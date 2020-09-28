%% Parse out valid file paths
function validFilePaths = ParseValidFilePath(origFilePaths)
    validFilePaths = cell(length(origFilePaths),1);
    index = 1;
    for i = 1:length(origFilePaths)
        if ~isempty(origFilePaths{i})
            validFilePaths{index} = origFilePaths{i};
            index = index + 1;
        end
    end
    if 1 == index
        validFilePaths = [];
    else
        validFilePaths = validFilePaths(1:index-1);
    end
end