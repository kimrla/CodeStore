function newindex =unidirectionsort ( qnrb,edgeindex,quadindex,startedgeindex,startquadindex,newindex,testnewindex )
% Sort a series of points in unidirection to a continuous curve.

% Input:
%   qnrb:Quad-nurbs structure.
%   edgeindex:Indices of intersection-edge corresponding to qnrb.qedges.
%   quadindex:Indices of quads corresponding to each edge in edgeindex
%       corresponding to qnrb.quad.It's a cell array.
%   startedgeindex:Index of the first intersection-edge corresponding to
%       qnrb.qedges,usually is edgeindex(1).
%   startquadindex:Index of the first quad which is connected with
%       startedgeindex,usually is quadindex{1}.
%   newindex:Index of the first intersection-edge corresponding to
%       edgeindex, usually is 1.
%   testnewindex:testnewindex=[testnewindex,newindex];testnewindex=unique(testnewindex);To
%       check if the goal has existed in newindex.

% Output:
%   newindex:Indices of intersection points after sorting them
%       corresponding to the original inputting intersections(edgeindex), which means
%       newpts=pts(newindex,:).

while (~isempty(startedgeindex))
        qd2edge=quad2edge(qnrb,startquadindex);
        % This means that there is only one intersection point in the edge being
        % considered and other intersections in this edge are ignored.
        tem=logical(qd2edge~=startedgeindex);
        qd2edge=qd2edge(tem);
        for i=1:3
            tem11=find(edgeindex==qd2edge(i));
            % If there are more than 2 intersection-edges in one quad
            % element and tem1 has been added into newindex, then tem1 
            % mustn't be considered for the second time.           
            if (~isempty(tem11))
                for kk=1:length(tem11)
                    tem1=tem11(kk);
                % Ensure each intersection-edge is only passed for one time.
                    tem8=sum(logical(testnewindex==tem1));
                    if (tem8==0)
                        secondedgeindex=qd2edge(i);
                        break;
                    end
                end
                if (tem8==0)
                    break;
                end
            end
        end
        if (~isempty(tem11) && tem8==0)
            newindex=[newindex,tem1];
            testnewindex=[testnewindex,newindex];
            testnewindex=unique(testnewindex);
            startedgeindex=secondedgeindex;
            tem=quadindex{tem1};
            if (length(tem)==2)
                tem_=tem~=startquadindex;
                startquadindex=tem(tem_);
            elseif (length(tem)==1)
                startedgeindex=[];
            end
        else
            startedgeindex=[];
        end 
        % If the intersection is a closed loop and the first point is not
        % coincident with the last one.(eg. if qnrb1 is the test
        % surface).And if so, the first element in newindex is NOT the same as
        % the last one.
        numuniq=length(unique(newindex));
        if (numuniq~=length(newindex))
            break;
        end
end



end

