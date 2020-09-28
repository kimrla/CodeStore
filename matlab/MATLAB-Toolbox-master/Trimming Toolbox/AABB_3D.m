function [pt1,pt2,center,id]=AABB_3D(mode,varargin)

% Construct 3D AABB box from scattered points or a NURBS/Bezier model,
% comparing with AABB_2D.m.

% Input:
%   mode: 'Sactter' for scattered points or 'Vertex' for NURBS/Bezier model
%   ptx,pty,ptz: The x,y,z coordinates of points, or the scattered points
%     pts ([numpnts,3]).
%   model: NURBS/bezier model structure.
% Output:
%   pt1: The left-down point of the box whose coordinate is (Xmin,Ymin,Zmin).
%   pt2: The right-up point of the box whose coordinate is (Xmax,Ymax,Zmax).
%   center: The centroid point of the box.
%   id: The 2 indices of pt1,pt2 corresponding to the scattered points,
%       or corresponding to the id of all the boundary box's vertices.

% If the input is a model, then construct its boundary box (polyhedron)
% firstly and then transform it to the AABB box using the vertices of the
% boundary box.

if strcmp(mode,'Scatter')
    if nargin==4
        pts=[varargin{1}(:),varargin{2}(:),varargin{3}(:)];
    elseif nargin==2
        pts=varargin{1};
    end
elseif strcmp(mode,'Vertex')
    model=varargin{1};   
    % Construct AABB of a NURBS/Bezier patch
    boundbox=nrbboundary(model);
    % In nrbboundary.m, if the model is NURBS, then firstly transform it to
    % several Bezier patch and the output is a cell structure.
    if strcmp(model.form,'B-NURBS')
        vertextract=@(model)(model.points(unique(model.id),:));
        tem1=cellfun(vertextract,boundbox,'UniformOutput',false);
        pts=cell2mat(tem1(:));
        pts=unique(pts,'rows');
    elseif strcmp(model.form,'Bezier')
        pts=boundbox.points(unique(boundbox.id),:);
        pts=unique(pts,'rows');
    end
else
    error('The mode should be ''Scatter'' or ''Vertex''');
end

% Construct AABB from scattered points    
[minpts,minid]=min(pts);
[maxpts,maxid]=max(pts);
pt1=minpts;
pt2=maxpts;
center=(pt1+pt2)/2;
if minid(1)==minid(2) && minid(2)==minid(3)
    id=minid(1);
else
    id=0;
end
if maxid(1)==maxid(2) && maxid(2)==maxid(3)
    id=maxid(1);
else
    id=0;
end    

end

%% demo
% pts=rand(10,3)*10;
% [pt1,pt2,center]=AABB_3D('Scatter',pts);
% % 8 vertices of the box.
% pt=[pt1(1) pt1(2) pt1(3)
%     pt2(1) pt1(2) pt1(3)
%     pt2(1) pt2(2) pt1(3)
%     pt1(1) pt2(2) pt1(3)
%     pt1(1) pt1(2) pt2(3)
%     pt2(1) pt1(2) pt2(3)
%     pt2(1) pt2(2) pt2(3)
%     pt1(1) pt2(2) pt2(3)];
% id=[1,2
%     2,3
%     3,4
%     4,1
%     5,6
%     6,7
%     7,8
%     8,5
%     1,5
%     2,6
%     3,7
%     4,8];
% figure;
% plot3(pts(:,1),pts(:,2),pts(:,3),'b*','MarkerSize',6);
% axis equal;hold on;    
% for i=1:length(id)   
%     plot3(pt(id(i,:),1),pt(id(i,:),2),pt(id(i,:),3),'LineWidth',2);
% end

%% demo
% srf=nrbtestsrf;
% box=nrbboundary(srf);
% figure;nrbctrlplot(srf)
% hold on;axis equal;
% k=3;
% trimesh(box{k}.id,box{k}.points(:,1),box{k}.points(:,2),box{k}.points(:,3));
% k=7;
% trimesh(box{k}.id,box{k}.points(:,1),box{k}.points(:,2),box{k}.points(:,3));
% [pt1,pt2,center]=AABB_3D('Vertex',srf);
% pt=[pt1(1) pt1(2) pt1(3)
%     pt2(1) pt1(2) pt1(3)
%     pt2(1) pt2(2) pt1(3)
%     pt1(1) pt2(2) pt1(3)
%     pt1(1) pt1(2) pt2(3)
%     pt2(1) pt1(2) pt2(3)
%     pt2(1) pt2(2) pt2(3)
%     pt1(1) pt2(2) pt2(3)];
% id=[1,2
%     2,3
%     3,4
%     4,1
%     5,6
%     6,7
%     7,8
%     8,5
%     1,5
%     2,6
%     3,7
%     4,8];  
% for i=1:length(id)   
%     plot3(pt(id(i,:),1),pt(id(i,:),2),pt(id(i,:),3),'LineWidth',2);
% end
% 
% 
% 

%% demo
% srf=nrbtestsrf;
% box=nrbboundary(srf);
% bzr=nrb2bzr(srf);
% figure;nrbctrlplot(srf)
% hold on;axis equal;
% k=3;
% trimesh(box{k}.id,box{k}.points(:,1),box{k}.points(:,2),box{k}.points(:,3));
% k=7;
% trimesh(box{k}.id,box{k}.points(:,1),box{k}.points(:,2),box{k}.points(:,3));
% [pt1,pt2,center]=AABB_3D('Vertex',bzr{7});
% pt=[pt1(1) pt1(2) pt1(3)
%     pt2(1) pt1(2) pt1(3)
%     pt2(1) pt2(2) pt1(3)
%     pt1(1) pt2(2) pt1(3)
%     pt1(1) pt1(2) pt2(3)
%     pt2(1) pt1(2) pt2(3)
%     pt2(1) pt2(2) pt2(3)
%     pt1(1) pt2(2) pt2(3)];
% id=[1,2
%     2,3
%     3,4
%     4,1
%     5,6
%     6,7
%     7,8
%     8,5
%     1,5
%     2,6
%     3,7
%     4,8];  
% for i=1:length(id)   
%     plot3(pt(id(i,:),1),pt(id(i,:),2),pt(id(i,:),3),'LineWidth',2);
% end
% 
% 
% 





