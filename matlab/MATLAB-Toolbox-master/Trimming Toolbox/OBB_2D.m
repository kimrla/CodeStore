function [ vert,center,sa,id ] = OBB_2D( varargin )
% Construct 2D OBB box from scattered points or the x/y coordinates. Use
% 'option' to change the method of constructing the box, which includes
% 'covariance'(default) for calculating the matrix or 'boundary' for using the
% MATLAB function boundary.m.

% Input:
%     ptx,pty: The x and y coordinates of points, or the scattered points pts.
%     option: 'covariance'(default) or 'boundary'.
% Output:
%     vert: The coordinates of the four vertices of box, which are
%       sorted.(size(vert)=[2,4])
%     center: The centroid point of the box,which is the intersection of
%       diagonals for both 2 kinds.
%     sa: The 2 separated axes of the box, represented by [2*2] matrix 
%       and each COLUMN vector stands one sa,which has been normalized.
%     id: The 4 indices of [lxmin,lxmax,lymin,lymax] corresponding to the input points.


if nargin==3
    pts=[varargin{1},varargin{2}];
    option=varargin{3};
elseif nargin==2
    if ischar(varargin{2})
        pts=varargin{1};
        option=varargin{2};
    else
        pts=[varargin{1},varargin{2}];
        option='covariance';
    end
else
    pts=varargin{1};
    option='covariance';
end

center=sum(pts)/length(pts);

if strcmp(option,'boundary')
    k=boundary(pts,0); % k: The first element is the same as the last one
    numk=length(k);
    s=10e6;
    for i=1:numk-1
        j=i+1;
        pt1=pts(k(i),:);
        pt2=pts(k(j),:);
        pt21=pt2-pt1;
        %theta: The rotation angle from the x-axis from global sys to local sys, which is represented by rad.
        theta=anglereverse(pt21);
        % Although cordtrans.m is for 2D, the input must be constructed as 3D coords.
        [trans,npt]=cordtrans([pts,zeros(length(pts),1)],[pt1,0],theta);
        [lymax,idymax]=max(npt(:,2));
        [lxmax,idxmax]=max(npt(:,1));
        [lxmin,idxmin]=min(npt(:,1));
        % Find the minimum area of the box.
        if s>(lxmax-lxmin)*lymax
            lxmin_=lxmin;
            lxmax_=lxmax;
            lymax_=lymax;
            trans_=trans;
            idxmin_=idxmin;
            idxmax_=idxmax;
            idymin_=k(i);
            idymax_=idymax;
            s=(lxmax-lxmin)*lymax;
        end
    end
    % Construct 3D homogeneous coords.
    lpts=[lxmin_,0,0,1;
        lxmax_,0,0,1;
        lxmax_,lymax_,0,1;
        lxmin_,lymax_,0,1];
    lpts=lpts';
    gpts=trans_\lpts;
    vert=gpts([1,2],:);
    sa(:,1)=vert(:,3)-vert(:,2);
    sa(:,2)=vert(:,2)-vert(:,1);
    sa=vecnorm(sa);
    id=[idxmin_,idxmax_,idymin_,idymax_];
    center=(vert(:,1)+vert(:,3))/2;
 
elseif strcmp(option,'covariance')
    covpts=cov(pts);
    [eigpts,~]=eig(covpts);
    sa=vecnorm(eigpts);
    if sum(cross([sa(:,1);0],[sa(:,2);0]))<0
        sa(:,2)=-sa(:,2);
    end
    theta=anglereverse(sa(:,1));
    [trans,npt]=cordtrans([pts,zeros(length(pts),1)],[center,0],theta);
    [lymax,idymax]=max(npt(:,2));
    [lymin,idymin]=min(npt(:,2));
    [lxmax,idxmax]=max(npt(:,1));
    [lxmin,idxmin]=min(npt(:,1));
    
    lpts=[lxmin,lymin,0,1;
        lxmax,lymin,0,1;
        lxmax,lymax,0,1;
        lxmin,lymax,0,1];
    lpts=lpts';
    gpts=trans\lpts;
    vert=gpts([1,2],:);
    id=[idxmin,idxmax,idymin,idymax];
    center=(vert(:,1)+vert(:,3))/2;
    
else
    warning('The method of constructing OBB_2D box is wrong.');
end


end


%% demo
% pt=rand(10,2)*10;
% [ vert,center,sa,id ] = OBB_2D( pt,'covariance' );vert(:,end+1)=vert(:,1);
% figure;hold on;
% plot(pt(:,1),pt(:,2),'b*','MarkerSize',6);
% plot(vert(1,:),vert(2,:),'r','LineWidth',2);
% quiver(center(1),center(2),sa(1,1)*10,sa(2,1)*10);
% quiver(center(1),center(2),sa(1,2)*10,sa(2,2)*10);
% 
% [ vert,center,sa,id ] = OBB_2D( pt,'boundary' );vert(:,end+1)=vert(:,1);
% plot(vert(1,:),vert(2,:),'g','LineWidth',2);
% quiver(center(1),center(2),sa(1,1)*10,sa(2,1)*10);
% quiver(center(1),center(2),sa(1,2)*10,sa(2,2)*10);
% axis equal;




