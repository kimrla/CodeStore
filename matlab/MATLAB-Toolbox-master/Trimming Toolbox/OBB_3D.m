function varargout=OBB_3D(varargin)

% Construct 3D OBB box from scattered points or the x/y coordinates. Use
% 'option' to change the method of constructing the box, which includes
% 'covariance'(default) of calculating the matrix or 'boundary' of
% transforming NURBS/Bezier boundary box to OBB.

% Input:
%     ptx,pty: The x and y coordinates of points, or the scattered points pts.
%     option: 'covariance'(default) or 'boundary'.
%     box: Control boundary box, a structure array, including 2D polygon
%       or 3D polyhedrons. The input box structure includes 
%       box.points/.id/.sa, see nrbboundary.m and ctrl2box.m.
%     ifmin:  If it is to find the minimum area/volumn box.The default is []
% Output:
%     vert: The coordinates of the 8 vertices of box, which are
%       sorted.(size(vert)=[8,3])
%     center: The centroid point of the box,which is the intersection of
%       diagonals for both 2 kinds.
%     sa: The 3 separated axes of the box, represented by [3*3] matrix 
%       and each COLUMN vector stands one sa,which has been normalized.
%     newbox: OBB box structure after calling ctrl2box.m

% Notcie that this function calls ctrl2box.m, so if mode=='boundary',the
% default OBB box is NOT the minimum volume. If want to find the minimum
% box, then find the function ctrl2box.m and use ifmin.
option=varargin{1};
if strcmp(option,'covariance')
    if nargin==2
        pts=varargin{2};
    elseif nargin==4
        pts=[varargin{2}(:),varargin{3}(:),varargin{4}(:)];        
    end
elseif strcmp(option,'boundary')
    if nargin==2
        box=varargin{2};
        ifmin=[];
    elseif nargin==3
        box=varargin{2};
        ifmin='min';
    end
else
    error('The first input should be ''covariance'' or ''boundary''');
end

% Only process one boundary box each time.
if strcmp(option,'boundary')
    if isempty(ifmin)
        newbox=ctrl2box(box);
    else
        newbox=ctrl2box(box,ifmin);
    end
    if nargout==3
        vert=newbox.points;
        vert=vert(:,1:3)./repmat(vert(:,4),1,3);
        sa=newbox.sa;
        sa=sa';% Each COLUMN means a sa.
        center=(vert(1,:)+vert(7,:))/2;
%         % 12 edges of the box
%         newbox.id=[1,2
%              2,3
%              3,4
%              4,1
%              1,5
%              2,6
%              3,7
%              4,8
%              5,6
%              6,7
%              7,8
%              8,5];              
    end
elseif strcmp(option,'covariance')    
    center=sum(pts)/length(pts);
    covpts=cov(pts);
    [eigpts,~]=eig(covpts);
    sa=vecnorm(eigpts);% Each COLUMN means a sa.
    if dot(cross(sa(:,1),sa(:,2)),sa(:,3))<0
        sa(:,3)=-sa(:,3);
    end
    mt=cordtrans3D(center,sa(:,1),sa(:,2),sa(:,3));
    tempts=pts';
    tempts(4,:)=1;% Homogeneous coords.
    npts=mt*tempts;
    npts=npts(1:3,:)./repmat(npts(4,:),3,1);
    npts=npts';% Non-homogeneous coords.
    lzmax=max(npts(:,3));
    lzmin=min(npts(:,3));    
    lymax=max(npts(:,2));
    lymin=min(npts(:,2));
    lxmax=max(npts(:,1));
    lxmin=min(npts(:,1));
    pt1=[lxmin,lymin,lzmin];
    pt2=[lxmax,lymax,lzmax];
    % 8 vertices of the box.
    pt=[pt1(1) pt1(2) pt1(3)
        pt2(1) pt1(2) pt1(3)
        pt2(1) pt2(2) pt1(3)
        pt1(1) pt2(2) pt1(3)
        pt1(1) pt1(2) pt2(3)
        pt2(1) pt1(2) pt2(3)
        pt2(1) pt2(2) pt2(3)
        pt1(1) pt2(2) pt2(3)];
    pt(:,4)=1;% Homogeneous coords.
    pt=pt';
    gpts=mt\pt;
    vert=gpts(1:3,:)./repmat(gpts(4,:),3,1);
    vert=vert';
    center=(vert(1,:)+vert(7,:))/2;

else
    warning('The method of constructing OBB_2D box is wrong.');
end

if nargout==1
    varargout{1}=newbox;
elseif nargout==3
    varargout{1}=vert;
    varargout{2}=center;
    varargout{3}=sa;
end

end

%% demo
% clear;
% pt=rand(10,3)*20;
% [ vert,center,sa] = OBB_3D( 'covariance',pt);
% figure;hold on;
% plot3(pt(:,1),pt(:,2),pt(:,3),'b*','MarkerSize',6);
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
%     plot3(vert(id(i,:),1),vert(id(i,:),2),vert(id(i,:),3),'LineWidth',2);
% end
% quiver3(center(1),center(2),center(3),sa(1,1)*10,sa(2,1)*10,sa(3,1)*10);
% quiver3(center(1),center(2),center(3),sa(1,2)*10,sa(2,2)*10,sa(3,2)*10);
% quiver3(center(1),center(2),center(3),sa(1,3)*10,sa(2,3)*10,sa(3,3)*10);
% axis equal;view(3);

%% demo
% clear;
% srf=nrbtestsrf;
% box=nrbboundary(srf);
% figure;nrbctrlplot(srf)
% hold on;
% k=4;
% trimesh(box{k}.id,box{k}.points(:,1),box{k}.points(:,2),box{k}.points(:,3));
% boundbox=box{k};
% [ vert,center,sa] = OBB_3D( 'boundary',boundbox,'min');
% % figure;hold on;
% plot3(vert(:,1),vert(:,2),vert(:,3),'b*','MarkerSize',6);
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
%     plot3(vert(id(i,:),1),vert(id(i,:),2),vert(id(i,:),3),'LineWidth',2);
% end
% quiver3(center(1),center(2),center(3),sa(1,1)*10,sa(2,1)*10,sa(3,1)*10,'LineWidth',2);
% quiver3(center(1),center(2),center(3),sa(1,2)*10,sa(2,2)*10,sa(3,2)*10,'LineWidth',2);
% quiver3(center(1),center(2),center(3),sa(1,3)*10,sa(2,3)*10,sa(3,3)*10,'LineWidth',2);
% axis equal;view(3);



