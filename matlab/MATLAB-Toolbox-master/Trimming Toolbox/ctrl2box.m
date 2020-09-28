function newbox=ctrl2box(boundbox,varargin)

% Transform the control boundary polyhedrons or polygons of a NURBS/Bezier
% model to OBB box. (2D or 3D). This function process one boundary box each
% time, instead the cell array of several boxes.

% Input:
%   boundbox: Control boundary box, a structure array, including 2D polygon
%       or 3D polyhedrons. The input box structure includes 
%       box.points/.id/.sa, see nrbboundary.m.
%   varargin{1}: If it is to find the minimum area/volumn box.
% Output:
%   newbox: New box structure of OBB after the transformation.
%       box.points/id/sa. Notice the box.points is homogeneous coords,
%       which means size(points)=[numpoints,4].

% For 2D, the box is minimum; For 3D, the box is NOT always minimum,
% depending on if there's the variation mode:'min'.

points=boundbox.points;
id=boundbox.id;
sa=boundbox.sa;% Number of ROWS equal the number of axes.

[num,dim]=size(id); % num is a GLOBAL variation
if dim==2 % 2D polygon->boundary box 
    pts=points(:,[1,2]);
    [vert,sa2d]=minobb2d(pts,id);
    newbox=boxconstruct_2D('OBB',vert,sa2d);
   
elseif dim==3 % 3D polyhedron-> boundary body
    if nargin==1          
        [vert3,sa2d,~]=obb3d(points,id,1,sa);       
         newbox.points=vert3'; 
         newbox.sa=[sa2d(1:3,:)';sa(1,:)];% Each ROW means a sa.
         % 12 edges of the box
         newbox.id=[1,2
             2,3
             3,4
             4,1
             1,5
             2,6
             3,7
             4,8
             5,6
             6,7
             7,8
             8,5];     
    elseif strcmp(varargin{1},'min')
        vol=10e8;
        for ii=1:num
            [vert3,sa2d,vol_]=obb3d(points,id,ii,sa);
            if vol_<vol
                vert3_=vert3;
                sa2d_=sa2d;
                ii_=ii;
            end
        end
        newbox.points=vert3_';
        newbox.sa=[sa2d_(1:3,:)';sa(ii_,:)];
        % 12 edges of the box
        newbox.id=[1,2
             2,3
             3,4
             4,1
             1,5
             2,6
             3,7
             4,8
             5,6
             6,7
             7,8
             8,5];                   
    else
        error('The mode is wrong');
    end

end
end

    function [vert,tsa]=minobb2d(pts,id)
        % Transform polygon to minimum 2D OBB
        s=10e6;
        [temnum,~]=size(id); % num is a GLOBAL variation
        for i=1:temnum 
            % size(pts)=[dim,2], 2D coords without the 3rd element, which means pts=points(:,[1,2])
            pt1=pts(id(i,1),:);
            pt2=pts(id(i,2),:);
            pt21=pt2-pt1;
            %theta: The rotation angle from the x-axis from global sys to local sys, which is represented by rad.
            theta=anglereverse(pt21);
            % Although cordtrans.m is for 2D, the input must be constructed as 3D coords.
            [trans,npt]=cordtrans([pts,zeros(length(pts),1)],[pt1,0],theta);
            [lymax,~]=max(npt(:,2));
            [lxmax,~]=max(npt(:,1));
            [lxmin,~]=min(npt(:,1));
            % Find the minimum area of the box.
            if s>(lxmax-lxmin)*lymax
                lxmin_=lxmin;
                lxmax_=lxmax;
                lymax_=lymax;
                trans_=trans;
                s=(lxmax-lxmin)*lymax;
            end
        end
        % Construct 4D homogeneous coords.
        lpts=[lxmin_,0,0,1;
            lxmax_,0,0,1;
            lxmax_,lymax_,0,1;
            lxmin_,lymax_,0,1];
        lpts=lpts'; % 4 vertices in local coord-sys.
        gpts=trans_\lpts; % gpts means homogeneous coords, comparing with vert.
        % size(vert)=[3,4]
        vert=gpts([1,2,3],:);% vert means 4 vertices of the box. In this place, vert means 3D coords.
        tsa(:,1)=vert(:,3)-vert(:,2);% size(tsa)=[3,2]
        tsa(:,2)=vert(:,2)-vert(:,1);
        tsa=vecnorm(tsa); 
    end

    function [vert3,sa2d,vol]=obb3d(points,id,i,temsa)
         % The ith patch and the 1st edge in this patch.
         temid=unique(id(:));
         boundpts=points(temid,:);% boundary vertices in global coord-sys (3D elements).
         pt1=points(id(i,1),:);
         pt2=points(id(i,2),:);
         x=pt2-pt1;    
         z=temsa(i,:);
         y=cross(z,x);
         mt=cordtrans3D(pt1,x,y,z);
         % Transform glocal coords to local coords.
         boundpts=boundpts';
         boundpts(4,:)=1;
         lboundpts=mt*boundpts;
         lboundpts=lboundpts';
         lboundpts(:,4)=[];
         zmax=max(lboundpts(:,3));
         zmin=min(lboundpts(:,3));
         tempts=lboundpts(:,[1,2]);
         kk=convhull(tempts);
         kk=[kk(1:end-1),kk(2:end)];
         [temvert,sa2d]=minobb2d(tempts,kk);
         % vert3: 8 vertices of the box. sa2d: 2 separated axes.Both are the homogeneous in
         % local coord-sys
         temvert(3,:)=zmin;
         temvert(4,:)=1;
         vert3=temvert;
         temvert(3,:)=zmax;
         vert3=[vert3,temvert];
         sa2d(4,:)=1;% size(sa2d)=[4,2]
%          center=(vert3(:,1)+vert3(:,7))/2;
         % Transform from local sys to global sys.
         vert3=mt\vert3;
         sa2d=mt\sa2d;
         center=mt\[0;0;0;1];
         sa2d=sa2d-repmat(center,1,2);
         % Calculate the volume of the box
         vol=norm(temvert(:,1)-temvert(:,2))*norm(temvert(:,3)-temvert(:,2))*(zmax-zmin);
    end
%%
% close all;
% srf=nrbtestsrf;
% box=nrbboundary(srf);
% figure;nrbctrlplot(srf)
% hold on;
% k=4;
% trimesh(box{k}.id,box{k}.points(:,1),box{k}.points(:,2),box{k}.points(:,3));
% %%
% boundbox=box{k};
% newbox=ctrl2box(boundbox);
% %%
% pts=newbox.points;
% id=newbox.id;
% for i=1:12
%     plot3(pts(id(i,:),1),pts(id(i,:),2),pts(id(i,:),3),'ko-','MarkerSize',4,'Linewidth',2);
% end
% %%
% newbox=ctrl2box(boundbox,'min');
% pts=newbox.points;
% id=newbox.id;
% for i=1:12
%     plot3(pts(id(i,:),1),pts(id(i,:),2),pts(id(i,:),3),'ro-','MarkerSize',4,'Linewidth',2);
% end
% 





