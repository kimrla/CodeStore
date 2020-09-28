function newindex2=indexmerge(pts,quadindex,newindex1)
% Merge parts of newindex1 into one series or two series of parts. This
% prodedure includes merging isolated points in newindex1 into their
% corresponding parts. 

% Input:
%   qnrb: Quad-nurbs structure.
%   pts: Parameter coordinates of the intersection points.
%   edgeindex: Indices of edges for each intersection points corresponding
%       to qnrb.qedges. 
%   quadindex: Indices of quads for each intersection edges corresponding
%       to qnrb.quad. 
%   newindex1: Indices of several curve parts in parameter domain after sorting
%       preliminarily, which includes all the points in pts/edgeptindex,see
%       newindex_{}.

%   Output:
%   newindex2: Indices of par-curve, which is a 1*1 cell array or 1*2 cell
%       array, and there should not be isolated points.

% A continuous intersection curve in physical domain may be separated into
% 2 parts in parameter domain (eg. revolving sphere surface), and
% furthermore, in parameter domain, a continuous par-intersection curve may
% be separated into several parts or isolated points because: 1. There are
% more than one intersection points in one edge and 2. There are more than
% 2 points in one quad element.

numnewindex1=length(newindex1);
isolatepoint=[];
nonisolatepoint={[]};
% Divide newindex1{i} into 2 kinds of groups: Isolated points group and
% non-isolated points group.
for i=1:numnewindex1
    numindex=length(newindex1{i});
    if (numindex==1)
        isolatepoint=[isolatepoint,newindex1{i}];% vector array
    else
        nonisolatepoint(end+1)=newindex1(i);% cell array
    end
end

nonisolatepoint(1)=[];
newindex2=nonisolatepoint;
numisolatepoint=length(isolatepoint);
numnonisolatepoint=length(nonisolatepoint);

% Insert all the isolated points into their corresponding curve parts.
if (numisolatepoint~=0)
    for i=1:numisolatepoint
        indexpt3=isolatepoint(i);
        pt3=pts(indexpt3,:);
        % Search the corresponding quads of the isolate point and the
        % existing line conncected by pt1 and pt2.
        quads=quadindex{indexpt3};
        for j=1:length(quads)
            k=quads(j);% k means the isolate point is in the kth quad.
            for jj=1:numnonisolatepoint% jj is the jjth non-isolated curve parts
                testpoint=nonisolatepoint{jj};
                for jjj=1:length(testpoint)-1% jjj is the jjjth element in one noisolatepoint
                    quad1=quadindex{testpoint(jjj)};% testpoint(jjj) is the testpoint(jjj)th intersection point.
                    quad2=quadindex{testpoint(jjj+1)};
                    goalquad=intersect(quad1,quad2);
                    if (goalquad==k)
                        pt1=pts(testpoint(jjj),:);
                        pt2=pts(testpoint(jjj+1),:);
                        break;
                    end
                end
                if (goalquad==k)
                    break;
                end
            end
            if (goalquad==k)
                break;
            end
        end
        [~,sortuv]=pointinsert(pt1,pt2,pt3);   
        if (isequal(sortuv,[1,2,3]))
            newindex2{jj}=[newindex2{jj}(1:jjj-1),newindex2{jj}(jjj),newindex2{jj}(jjj+1),indexpt3,newindex2{jj}(jjj+2:end)];       
        elseif (isequal(sortuv,[1,3,2]))
            newindex2{jj}=[newindex2{jj}(1:jjj-1),newindex2{jj}(jjj),indexpt3,newindex2{jj}(jjj+1),newindex2{jj}(jjj+2:end)];           
        else
            newindex2{jj}=[newindex2{jj}(1:jjj-1),indexpt3,newindex2{jj}(jjj),newindex2{jj}(jjj+1),newindex2{jj}(jjj+2:end)];
        end
    end
end
    

% If there are more than 2 cell arrays in nonisolatepoint{}, all the cell
% arrays need to be merged into one or two cells. And if there are 2 cells,
% it needs to be checked if the 2 cells can be merged into one one cell or
% a closed loop.
if (numnonisolatepoint>=2)
    for i=1:numnonisolatepoint-1
        temindex1(1)=newindex2{i}(1);
        temindex1(2)=newindex2{i}(end);
        temindex2(1)=newindex2{i+1}(1);
        temindex2(2)=newindex2{i+1}(end);
        for j=1:2
            index1=temindex1(j);
            for k=1:2
                index2=temindex2(k);
                quad1=quadindex{index1};
                quad2=quadindex{index2};
                quad3=unique([quad1,quad2]);
                if (length(quad3)==length(quad1)+length(quad2)-1)
                    if (j==1 && k==1)
                        newindex2{i}=[newindex2{i}(end:-1:1),newindex2{i+1}];    
                    elseif (j==1 && k==2)
                        newindex2{i}=[newindex2{i}(end:-1:1),newindex2{i+1}(end:-1:1)]; 
                    elseif (j==2 && k==1)
                        newindex2{i}=[newindex2{i},newindex2{i+1}]; 
                    else
                        newindex2{i}=[newindex2{i},newindex2{i+1}(end:-1:1)]; 
                    end              
                    newindex2(i+1)=[];
                    break;
                end
            end
            if (length(quad3)==length(quad1)+length(quad2)-1)
                break;
            end
        end
    end   
else                                            
    % Check if the curve is a closed loop.
    index1=newindex2{1}(1);
    index2=newindex2{1}(end);
    quad1=quadindex{index1};
    quad2=quadindex{index2};
    quad3=unique([quad1,quad2]);
    if (length(quad3)==length(quad1)+length(quad2)-1)
        newindex2{1}=[newindex2{1},newindex2{1}(1)];
    end    
end
        
end