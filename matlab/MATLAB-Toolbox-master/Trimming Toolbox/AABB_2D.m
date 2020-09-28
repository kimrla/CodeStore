function [ pt1,pt2,center,id ] = AABB_2D( varargin )
% Construct 2D AABB box from scattered points or the coordinates.

% Input:
%     ptx,pty: The x and y coordinates of points, or the scattered points pts ([numpnts,2]).
% Output:
%     pt1: The left-down point of the box whose coordinate is (Xmin,Ymin).
%     pt2: The right-up point of the box whose coordinate is (Xmax,Ymax).
%     center: The centroid point of the box.
%     id: The 2 indices of pt1 and pt2 corresponding to the input points.

if nargin==2
    pts=[varargin{1},varargin{2}];
else
    pts=varargin{1};
end

[minpts,minid]=min(pts);
[maxpts,maxid]=max(pts);

pt1=minpts;
pt2=maxpts;
center=(pt1+pt2)/2;

if minid(1)==minid(2)
    id=minid(1);
else
    id=0;
end
if maxid(1)==maxid(2)
    id=maxid(1);
else
    id=0;
end


end

%% demo
% pts=rand(10,2)*10;
% figure;plot(pts(:,1),pts(:,2),'b*','MarkerSize',6);axis equal;hold on;
% [pt1,pt2,center]=AABB_2D(pts);
% pt=[pt1;pt2(1) pt1(2);pt2;pt1(1) pt2(2);pt1];
% plot(pt(:,1),pt(:,2),'LineWidth',2);


