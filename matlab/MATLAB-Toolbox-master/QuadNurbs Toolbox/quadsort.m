function [newindex_,newpts,newpnts]=quadsort(qnrb,pts,pnts,edgeptindex)
% Sort discrete intersection points as a continuous curve in physical
% domain, while in parameter domain, the ONE continuous
% physical-intersection curve may be separated into 2 parts in degenerate
% surface, eg. Revolved Sphere surface.

% Input:
%   qnrb:Quad-nurbs surface structure.
%   pts:Parameter coordinates of the intersection points.
%   pnts:Physical coordinates of the intersection points.
%   edgeptindex:Edges corresponding to each intersection point, which are 
%       represented by indices of 2 ending points of the edge. See ed1 in [ed1, tri2, pts1, pts2, pnts1,
%       pnts2]=qnrbinter(qnrb1, qnrb2).
% Output:
%   newpts:New parameter coordinates of intersection points after sorting
%       them, which is a cell array corresponding to newindex_.
%   newpnts:New physical coordinates of intersection points after sorting
%       them,, which is a cell array corresponding to newindex_.
%   newindex_:Indices of intersection points after sorting them
%       corresponding to the original inputting intersections, which means
%       newpts=pts(newindex,:).

% If the curve in parameter domain is separated into 2 parts(eg.The
% revolved sphere surface's revolving-coincidence edge), then the output
% newindex is a 1*2 cell array, otherwise it's a 1*1 cell array.

edges=qnrb.qedges;
ed2qd=qnrb.ed2qd;
edges=sort(edges,2);
edgeptindex=sort(edgeptindex,2);
testnewindex=[];

% Remove dumplicate intersections (dumplicate points means the same points
% in both parameter and physical domain,but those points which belongs to the same intersection-edge are not included.).
[~,u,v]=unique(edgeptindex,'rows','stable');
num1=length(edgeptindex);
num2=length(u);
dumplicateedgeindex(1)={[]};
if(num1~=num2)
    removeindex=[];
    max1=max(v);
    for i=1:max1
        tem5=find(v==i);
        if (length(tem5)>1)
            dumplicateedgeindex(end+1)={tem5};
        end
    end
    dumplicateedgeindex(1)=[];
    % dumplicateedgeindex is a cell array and stores all the same edges' indices
    % corresponding to edgeptindex.
    numtem6=length(dumplicateedgeindex);
    for j=1:numtem6
        testindex=dumplicateedgeindex{j};
        numtestindex=length(testindex);
        for k=2:numtestindex
            tem7=norm(pts(testindex(1),:)-pts(testindex(k),:));
            if (tem7<1e-6)
                removeindex=[removeindex,testindex(k)];
            end
        end
    end
    if (~isempty(removeindex))
        edgeptindex(removeindex,:)=[];
        pts(removeindex,:)=[];
        pnts(removeindex,:)=[];
    end
end
       
n=length(edgeptindex);
logedge=logical(false(length(edges),1));
newindex=1;
newindex_{1}=[];
testnewindex=[testnewindex,newindex];
testnewindex=unique(testnewindex);

% Find the relationship between edgeindex,qnrb.qedges and qnrb.ed2qd.
for i=1:n
    logedge= find(edges(:,1)==edgeptindex(i,1) & edges(:,2)==edgeptindex(i,2));
    edgeindex(i)=logedge;
end
quadindex=ed2qd(edgeindex);

% Sorting intersections from the first point/edgeindex/dealquad.
startquadindex=quadindex{1};
startedgeindex=edgeindex(1);
m=length(startquadindex);

while (true)
    if (m==1)
    % The edge is connected with one quad.
        newindex =unidirectionsort ( qnrb,edgeindex,quadindex,startedgeindex,startquadindex,newindex,testnewindex );    
    else
    % The edge is connected with 2 quads, which means m=2.Sort the
    % intersections along with 2 directions and then merge them.
        for i=1:2
            temindex{i}=unidirectionsort(qnrb,edgeindex,quadindex,startedgeindex,startquadindex(i),newindex,testnewindex );
        end
        temindex2=temindex{1}(end:-1:1);
        % If the intersection curve is a closed loop and the first point is not
        % coincident with the last one, then [temindex2(end),temindex2] is the same as
        % [temindex{2},temindex{2}(1)] and the first element in it is NOT the same as the last
        % one. If there's only an isolated point in temindex, then
        % temindex{1}=temindex{2}=isolated point index.
        if (isequal([temindex2(end),temindex2],[temindex{2},temindex{2}(1)]))
            newindex=temindex2;
        else
            temindex2(end)=[];
            newindex=[temindex2,temindex{2}];
        end

    end

    newindex_{end+1}=newindex;
    testnewindex=[testnewindex,newindex];
    testnewindex=unique(testnewindex);
    numtestnewindex=length(testnewindex);

    % If there are more than one parts of the entire curve in the parameter
    % domain.Ohter parts may include only one isolated point.
    if (numtestnewindex~=n)
    % n=length(edgeindex) is the number of intersections after removing dumplicate points, 
    % while num1=length(edgeptindex) is the number of inputting intersections 
        for jj=1:numtestnewindex-1
            if ((testnewindex(jj)+1)~=testnewindex(jj+1))
                newindex=testnewindex(jj)+1;
                startedgeindex=edgeindex(newindex);
                startquadindex=quadindex{newindex};
                break;
            end
        end
        if (jj==numtestnewindex-1 && (testnewindex(jj)+1)==testnewindex(jj+1) )
            newindex=testnewindex(end)+1;
            startedgeindex=edgeindex(newindex);
            startquadindex=quadindex{newindex};
        end
        m=length(startquadindex);
    else
        break;
    end
end
    
newindex_(1)=[];    
% Merge all the parts in newindex_{} as soon as possible.
newindex_=indexmerge(pts,quadindex,newindex_);
    
if (nargout>1)
    for i=1:length(newindex_)
        newpts{i}=pts(newindex_{i},:);
        newpnts{i}=pnts(newindex_{i},:);    
    end
end

end
    
%% demo
% h0=1.8;
% 
% % Create a nurbs sphere
% L=[7,10,0; 5,1,10];
% center=[5,5,4];
% circ=nrbcirc(4, center, 0, pi);
% srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% 
% % Transform a nurbs surface into quadrangular representation
% qnrb1=nrb2quad(srf1, h0);
% qnrb2=nrb2quad(srf2, h0);
% 
% % Get approximated intersections by quadrangles and triangles
% [ed1, tri2, pts1, pts2, pnts1, pnts2]=qnrbinter(qnrb1, qnrb2);
% 
% % Plot the surface and related quadrangle and line
% figure; hold on;
% trisurf(qnrb1.quad, qnrb1.points(:,1), qnrb1.points(:,2), qnrb1.points(:,3));
% trisurf(qnrb2.tri, qnrb2.points(:,1), qnrb2.points(:,2), qnrb2.points(:,3));
% plot3(pnts1(:,1), pnts1(:,2), pnts1(:,3), 'r*');
% plot3(pnts2(:,1), pnts2(:,2), pnts2(:,3), 'ro');
% view(3); axis equal;
% 
% % Plot the parameric domain of surface 1
% figure; hold on;
% quadplot(qnrb1, 'r');
% plot(pts1(:,1), pts1(:,2), 'r*');
% axis equal;
% 
% % Plot the parameric domain of surface 2
% figure; hold on;
% triplot(qnrb2.tri, qnrb2.nodes(:,1), qnrb2.nodes(:,2), 'r');
% plot(pts2(:,1), pts2(:,2), 'r*');
% axis equal;
% 
% [newindex_,newpts,newpnts]=quadsort(qnrb1,pts1,pnts1,ed1);
% plot(newpts{1}(:,1),newpts{1}(:,2),'o-.','MarkerSize',2);
% 
% 
% 

    
    

    
    


