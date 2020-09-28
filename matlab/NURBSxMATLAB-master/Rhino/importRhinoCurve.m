%% function
% import Rhino 'List'
% 
% Xu Yi, 2019

%%
function [degree, knotVector, cvPt] = importRhinoCurve(filename)
fileID = fopen(filename);
stringTemp = '0'; % temp = fscanf(fid,'%s',1);
while ~strcmp( stringTemp, 'order' ) % Compare strings
    stringTemp = fscanf(fileID,'%s',1); % A = fscanf(fileID,formatSpec,sizeA)
end
fscanf(fileID,'%s',1); % ����'='
stringTemp = fscanf(fileID,'%s',1);
order = str2double(stringTemp);
degree = order -1;

while ~strcmp( stringTemp, 'cv_count' )
    stringTemp = fscanf(fileID,'%s',1); % A = fscanf(fileID,formatSpec,sizeA)
end
fscanf(fileID,'%s',1); % ����'='
stringTemp = fscanf(fileID,'%s',1);
cv_count = str2double(stringTemp);

while ~strcmp( stringTemp, 'Knot' )
    stringTemp = fscanf(fileID,'%s',1); % A = fscanf(fileID,formatSpec,sizeA)
end
fscanf(fileID,'%s',1); % ����'Vector'
fscanf(fileID,'%s',1); % ����'('
stringTemp = fscanf(fileID,'%s',1);
knotNum = str2double(stringTemp);
knotNum = knotNum+2; % Rhino�Ľڵ�������������(��һ�������һ��)
fgetl(fileID); % ��������
fgetl(fileID);

knotVector = zeros(1,knotNum);
i = 2;
fscanf(fileID,'%s',1); % ����'0'
while ~strcmp( stringTemp, 'Control' )
    stringTemp = fscanf(fileID,'%s',1);
    knotTemp = str2double(stringTemp);
    knotVector(i) = knotTemp;
    i = i+1;
    stringTemp = fscanf(fileID,'%s',1);
    knotMult = str2double(stringTemp);
    fgetl(fileID);
    if knotMult ~= 1
        for j = 1:knotMult-1
            knotVector(i) = knotTemp;
            i = i+1;
        end
    end
    stringTemp = fscanf(fileID,'%s',1); % ����index
end
knotVector(1) = knotVector(2);
knotVector(end) = knotVector(end-1);

fgetl(fileID); % ���� Control Points ����
fgetl(fileID);
cvPt = zeros(cv_count,3);
for i = 1:cv_count
    fscanf(fileID,'%s',1);
    if i <= 10
        fscanf(fileID,'%s',1);
    end
    for j = 1:3
        stringTemp = fscanf(fileID,'%s',1);
        stringTemp = strrep(stringTemp, '(', ''); % modifiedStr = strrep(origStr, oldSubstr, newSubstr) �ÿհ��滻'('����ȥ����'('��
        cvPtTemp = textscan(stringTemp,'%f');
        cvPt(i,j) = cell2mat(cvPtTemp);
    end
end

fclose(fileID);
end
