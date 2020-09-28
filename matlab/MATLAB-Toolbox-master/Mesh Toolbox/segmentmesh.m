function trimodel=segmentmesh(modeltree)

% Use the tree structure and the nodes obtained by bzrsegment('adaptive')
% to generate triangulr meshes.

% Input:
%   modeltree: Tree structure obtained by 'adaptive' bzrsegment.m
% Output:
%   trimodel: Structure array containing triangularion and
%       parametric/physical coordinates information.

% trimodel structure:
% trimodel.model:NURBS/Bezier model in the root node,the original input model.
% trimodel.knots: Par-coords of the original surface.
% trimodel.points: Phy-coords corresponding to par-coords.
% trimodel.triangulation: Triangulation topology.

idleaf=findmodeltree(modeltree,0,'child'); % Index of leaf node.
idroot=findmodeltree(modeltree,1,'depth'); % Index of root node.
% Par-coords and phy-coords.
knots=[];
pnts=[];
tri=[];

if ~isempty(idroot)
    trimodel.model=modeltree{idroot}.model;
else
    idroot=findmodeltree(modeltree,2,'depth');
    trimodel.model=modeltree(idroot);
end
% Control the decimal digits: roundn.Only available for surface instead of curve.
numleaf=length(idleaf);
for i=1:numleaf
    model=modeltree{idleaf(i)}.model;
    % Par-coords, every 4 points represent one patch's 4 vertices.
    [u,v]=meshgrid(unique(model.knots{1}),unique(model.knots{2}));
    u=u';v=v';
    knots=[knots;[u(:),v(:)]];
    % Phy-coords, each column represents one point
    pnts=[pnts,model.coefs(:,1,1),model.coefs(:,end,1),model.coefs(:,1,end),model.coefs(:,end,end)];
end
% Process the phy-coords.
pnts=pnts';% Each row represents one point.
if size(pnts,2)==4 % Homogeneous coords
    pnts=pnts./repmat(pnts(:,4),1,4); % Non-homogeneous coords.
    pnts(:,4)=[];
end
% Assemble all the par-knots and the corresponding phy-points.
% Points are sorted in each v-direction iso-par-line
knots=[knots(:,2),knots(:,1)]; % knot:[v,u]
[knots_,ia,~]=unique(knots,'rows'); % knot_:[v,u]
pnts_=pnts(ia,:);
% Points are sorted in each u-direction iso-par-line
temknots=[knots_(:,2),knots_(:,1)]; % temknot:[u,v]
[~,temia,temic]=unique(temknots,'rows');
% Process each patch: leaf node.
knots=[knots(:,2),knots(:,1)];
knots_=temknots;
for i=0:numleaf-1
    idmid1=[];idmid2=[];idmid3=[];idmid4=[];
    temidmid_=[];
    idu1=find(knots_(:,1)==knots(i*4+1,1) & knots_(:,2)==knots(i*4+1,2));
    idu2=find(knots_(:,1)==knots(i*4+2,1) & knots_(:,2)==knots(i*4+2,2));
    idu3=find(knots_(:,1)==knots(i*4+3,1) & knots_(:,2)==knots(i*4+3,2));
    idu4=find(knots_(:,1)==knots(i*4+4,1) & knots_(:,2)==knots(i*4+4,2));
    % Find if there are middle points at the u-direction edge of each patch.
    if idu2~=idu1+1 
        idmid1=[idmid1,idu1+1:idu2-1];
    end
    if idu4~=idu3+1 
        idmid2=[idmid2,idu3+1:idu4-1];
    end
    % Find if there are middle points at the v-direction edge of each patch.
    idu1_=temic(idu1);
    idu2_=temic(idu2);
    idu3_=temic(idu3);
    idu4_=temic(idu4);
    if idu3_~=idu1_+1
        temidmid_=[temidmid_,idu1_+1:idu3_-1];
        temidmid=temia(temidmid_);
        idmid3=[idmid3,temidmid'];
    end
    temidmid_=[];
    if idu4_~=idu2_+1
        temidmid_=[temidmid_,idu2_+1:idu4_-1];
        temidmid=temia(temidmid_);
        idmid4=[idmid4,temidmid'];
    end
    % If there are middle points, then insert new point using nrbeval.
    idmid=[idmid1,idmid2,idmid3,idmid4];
    if ~isempty(idmid)
        % Calculate the barycenter of all the boundary points in one patch
        %考虑用边界点的物理坐标计算形心，而不是用参数坐标
        boundpts=knots_([idmid,idu1,idu2,idu3,idu4],:);
        knotsnew=sum(boundpts)/length(boundpts);
        knots_=[knots_;knotsnew];
        idend=length(knots_);
        pntsnew=nrbeval(modeltree{idleaf(i+1)}.model,knotsnew);
        pnts_=[pnts_;pntsnew'];
        % Connect all the points in the 4 boundary of the patch with the new
        % inserting middle point.
        if isempty(idmid1)
            tri=[tri;idu1,idu2,idend];
        else
            temidmid=[idu1,idmid1,idu2];
            for j=1:length(temidmid)-1
                tri=[tri;temidmid(j),temidmid(j+1),idend];
            end
        end
        if isempty(idmid2)
            tri=[tri;idu3,idu4,idend];
        else
            temidmid=[idu3,idmid2,idu4];
            for j=1:length(temidmid)-1
                tri=[tri;temidmid(j),temidmid(j+1),idend];
            end
        end       
        if isempty(idmid3)
            tri=[tri;idu1,idu3,idend];
        else
            temidmid=[idu1,idmid3,idu3];
            for j=1:length(temidmid)-1
                tri=[tri;temidmid(j),temidmid(j+1),idend];
            end
        end        
        if isempty(idmid4)
            tri=[tri;idu2,idu4,idend];
        else
            temidmid=[idu2,idmid4,idu4];
            for j=1:length(temidmid)-1
                tri=[tri;temidmid(j),temidmid(j+1),idend];
            end
        end        
    else
        tri=[tri;idu1,idu2,idu3;idu2,idu3,idu4];
    end
        
end        

trimodel.knots=knots_;
trimodel.points=pnts_;
trimodel.triangulation=tri;

end

%%
% clear global;
% global NUMNODE DEPTH modeltree
% tol=0.5;
% srf=nrbtestsrf;
% bzr.form='Bezier';bzrsegment('Adaptive',srf,tol,1);
% % bzrsegment(bzr,tol,'Adaptive',1,0,1);
% % bzrsegment(bzrsrf,tol,'Adaptive',1);
% % bzrsegment(srf,tol,'Adaptive',1);
% trimodel=segmentmesh(modeltree);
% figure;
% trisurf(trimodel.triangulation,trimodel.points(:,1),trimodel.points(:,2),...
%     trimodel.points(:,3));
% axis equal;
% figure;
% hold on;
% knots=trimodel.knots;
% plot(knots(:,1),knots(:,2),'ro');
% triplot(trimodel.triangulation,knots(:,1),knots(:,2));
% axis equal;




