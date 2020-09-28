function logic=colldetect_3D(box1,box2)

% Collision detection of 2 boxes based on the projection theorem (for AABB 
% and OBB) or using the intersection detection among all edges and polygons
% (for NURBS control polygons).
% This function is to determine whether 2 polygons in 2D are intersected.

% Input:
%     box1,box2: Two boxes to be detected, which are represented by the
%       structure obtained by the function boxconstruct.
% Output:
%     logic: The logical scalar meaning 1: if the two boxes are intersected
%       or 0: if the two boxes are NOT intersected.

form1=box1.form;
form2=box2.form;

if strcmp(form1,'Sphere') || strcmp(form2,'Sphere') % If circle box exists
    if strcmp(form1,form2) % Sphere-Sphere
        center1=box1.center;
        center2=box2.center;
        r1=box1.radius;
        r2=box2.radius;
        if norm(center1-center2)>r1+r2
            logic=false;
        else
            logic=true;
        end
    else % Sphere-Polyhedrons( including rectangles)
        if strcmp(form1,'Sphere')
            center=box1.center;
            r=box1.radius;
            box=box2;
        else
            center=box2.center;% Row vector
            r=box2.radius;
            box=box1;
        end
        % All vertices of the boundary box and the 1st element is NOT the
        % same as the last one. size(pts)=[numpoints,3].
        pts=box.points;
        logic=false;
        dist=pts-repmat(center,length(pts),1);
        dist=sqrt(sum(dist.*dist,2));
        if min(dist)<=r
            % There's at leat 1 vertex inside the sphere.
            logic=true;
        else
            idface=box.idface;
            idface=idface(:,1:3);
            % Check each face of the polyhedron or rectangle.
            for i=1:length(idface)
                [~,~,dist]=interpoint2tri([pts(idface(i,1),:);...
                    pts(idface(i,2),:);pts(idface(i,3),:)],center);
                if dist<=r
                    logic=true;
                    break;
                end
            end          
        end    
    end
    
elseif  strcmp(form1,'POLYHEDRON') || strcmp(form2,'POLYHEDRON')% Polyhedrons exist
    % In colldetect_2D.m, the method is based on inpolygon to determine if
    % the vertices in one polygon is inside the other one, and if not, use
    % the separated-axis-projection method to determine furthermore.
    % In this function, considering there are many faces of a POLYHEDRON,
    % the separated-axis-method is unavailable. so the method is to use all
    % the edges of one POLYHEDRON to determine the intersection of each
    % face of the other POLYHEDRON, based on interpoint2tri.m.
    pts1=box1.points;
    pts2=box2.points;
    idedge=box1.idedge;
    idface=box2.face;
    logic=false;
    for i=1:length(pts1)
        for j=1:length(pts2)
            [~,~,par,~,inter]=interline2tri(pts2(idface(j,:),:),...
                pts1(idedge(i,:),:));
            if inter~=0 && (sum(par>=0 && par<=1)==3)
                logic=true;
                break;
            end
        end
        if inter~=0 && (sum(par>=0 && par<=1)==3)
            logic=true;
            break;
        end
    end
 
elseif  strcmp(form1,'OBB') || strcmp(form2,'OBB')% OBB exists and Polygons don't exist
    % Construct 15 separated axes.
    sa1=box1.sa;% Number of ROWS equal number of sa, obtained by boxconstruct_3D.m   
    sa2=box2.sa;
    sa1=sa1';sa2=sa2';% Number of COLUMNS equal number of sa
    sa3=cross(sa1,sa2);
    sa4=cross([sa1(:,[2,3]),sa1(:,1)],sa2);
    sa5=cross([sa1(:,3),sa1(:,[1,2])],sa2);
    sa=[sa1,sa2,sa3,sa4,sa5];
    sa=sa';% Number of ROWS equal number of sa
    % Remove repeated sa and sa=[0,0,0];
    sa=unique(sa,'rows');
    sa((sa(:,1)==0 & sa(:,2)==0 & sa(:,3)==0),:)=[];
    % All points of the two boxes' vertices   
    pts1=box1.points;
    pts2=box2.points;
    center1=(pts1(1,:)+pts1(7,:))/2;
    center2=(pts2(1,:)+pts2(7,:))/2;    
    center=center2-center1;
    % Length, width, height of the 2 boxes.    
    l1=(pts1(2,:)-pts1(1,:))/2;
    w1=(pts1(3,:)-pts1(2,:))/2;
    h1=(pts1(5,:)-pts1(1,:))/2;
    l2=(pts2(2,:)-pts2(1,:))/2;
    w2=(pts2(3,:)-pts2(2,:))/2;
    h2=(pts2(5,:)-pts2(1,:))/2;
    logic=true;    
    for i=1:size(sa,1)
        det1=abs(dot(center,sa(i,:))); 
        det2=max([abs(dot((l1+w1+h1),sa(i,:))),abs(dot((l1+w1-h1),sa(i,:))),...
            abs(dot((-l1+w1+h1),sa(i,:))),abs(dot((-l1+w1-h1),sa(i,:)))])+...
            max([abs(dot((l2+w2+h2),sa(i,:))),abs(dot((l2+w2-h2),sa(i,:))),...
            abs(dot((-l2+w2+h2),sa(i,:))),abs(dot((-l2+w2-h2),sa(i,:)))]);
        if det1>det2
            logic=false;
            break;
        end
    end
 
else % AABB-AABB
    logic=true;
    pts1=box1.points;
    pts2=box2.points;
    xspan1=[pts1(1,1),pts1(3,1)];
    xspan2=[pts2(1,1),pts2(3,1)];
    yspan1=[pts1(1,2),pts1(3,2)];
    yspan2=[pts2(1,2),pts2(3,2)];
    zspan1=[pts2(1,3),pts2(5,3)];
    zspan2=[pts2(1,3),pts2(5,3)];
    % At least one projection indicates NOT intersecting
    if xspan1(2)<xspan2(1) || xspan2(2)<xspan1(1) 
        logic=false;
    elseif yspan1(2)<yspan2(1) || yspan2(2)<yspan1(1) 
        logic=false;
    elseif zspan1(2)<zspan2(1) || zspan2(2)<zspan1(1)
        logic=false;
    end
            
end

%% demo
% clear;close all;
% pt1=rand(10,3)*20-10;
% pt2=rand(10,3)*20+10;
% % pt1=rand(10,3)*20;
% % pt2=rand(10,3)*20;
% [ vert,~,sa] = OBB_3D( 'covariance',pt1);
% box1=boxconstruct_3D('OBB',vert,sa);
% vert=box1.points;
% id=box1.idedge;
% center=(vert(1,:)+vert(7,:))/2;
% sa=box1.sa;
% %%
% figure;hold on;
% plot3(pt1(:,1),pt1(:,2),pt1(:,3),'b*','MarkerSize',6);
% for i=1:length(id)   
%     plot3(vert(id(i,:),1),vert(id(i,:),2),vert(id(i,:),3),'LineWidth',2);
% end
% quiver3(center(1),center(2),center(3),sa(1,1)*10,sa(1,2)*10,sa(1,3)*10);
% quiver3(center(1),center(2),center(3),sa(2,1)*10,sa(2,2)*10,sa(2,3)*10);
% quiver3(center(1),center(2),center(3),sa(3,1)*10,sa(3,2)*10,sa(3,3)*10);
% axis equal;view(3);
% %%
% [ vert_,~,sa_] = OBB_3D( 'covariance',pt2);
% box2=boxconstruct_3D('OBB',vert_,sa_);
% 
% log=colldetect_3D(box1,box2);
% 
% vert=box2.points;
% id=box2.idedge;
% center=(vert(1,:)+vert(7,:))/2;
% sa=box2.sa;
% %%
% plot3(pt2(:,1),pt2(:,2),pt2(:,3),'b*','MarkerSize',6);
% for i=1:length(id)   
%     plot3(vert(id(i,:),1),vert(id(i,:),2),vert(id(i,:),3),'LineWidth',2);
% end
% quiver3(center(1),center(2),center(3),sa(1,1)*10,sa(1,2)*10,sa(1,3)*10);
% quiver3(center(1),center(2),center(3),sa(2,1)*10,sa(2,2)*10,sa(2,3)*10);
% quiver3(center(1),center(2),center(3),sa(3,1)*10,sa(3,2)*10,sa(3,3)*10);
% axis equal;view(3);

%% demo
% clear;close all;
% figure;hold on;axis equal;view(3);
% srf=nrbtestsrf;
% box=nrbboundary(srf);
% k=4;
% boundbox=box{k};
% [ vert,~,sa] = OBB_3D( 'boundary',boundbox,'min');
% box1=boxconstruct_3D('OBB',vert,sa);
% vert=box1.points;
% id=box1.idedge;
% center=(vert(1,:)+vert(7,:))/2;
% sa=box1.sa;
% %%
% nrbplot(srf,[50,50],'light','on');
% trimesh(boundbox.id,boundbox.points(:,1),boundbox.points(:,2),...
%     boundbox.points(:,3));
% for i=1:length(id)   
%     plot3(vert(id(i,:),1),vert(id(i,:),2),vert(id(i,:),3),'LineWidth',2);
% end
% quiver3(center(1),center(2),center(3),sa(1,1)*10,sa(1,2)*10,sa(1,3)*10);
% quiver3(center(1),center(2),center(3),sa(2,1)*10,sa(2,2)*10,sa(2,3)*10);
% quiver3(center(1),center(2),center(3),sa(3,1)*10,sa(3,2)*10,sa(3,3)*10);
% 
% %%
% k=7;
% boundbox=box{k};
% [ vert_,~,sa_] = OBB_3D( 'boundary',boundbox,'min');
% box2=boxconstruct_3D('OBB',vert_,sa_);
% 
% log=colldetect_3D(box1,box2);
% 
% vert=box2.points;
% id=box2.idedge;
% center=(vert(1,:)+vert(7,:))/2;
% sa=box2.sa;
% %%
% for i=1:length(id)   
%     plot3(vert(id(i,:),1),vert(id(i,:),2),vert(id(i,:),3),'LineWidth',2);
% end
% quiver3(center(1),center(2),center(3),sa(1,1)*10,sa(1,2)*10,sa(1,3)*10);
% quiver3(center(1),center(2),center(3),sa(2,1)*10,sa(2,2)*10,sa(2,3)*10);
% quiver3(center(1),center(2),center(3),sa(3,1)*10,sa(3,2)*10,sa(3,3)*10);
% 
% 




