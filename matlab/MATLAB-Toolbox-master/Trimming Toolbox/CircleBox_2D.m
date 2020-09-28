function [ center,r ] = CircleBox_2D( varargin )
% Construct 2D Circle box from scattered points or the coordinates.

% Input:
%     ptx,pty: The x and y coordinates of points, or the scattered points pts.
% Output:
%     center: The centroid point of the circle box.
%     r: The radius of the circle.


if nargin==2
    pts=[varargin{1},varargin{2}];
else
    pts=varargin{1};
end

center=sum(pts)/length(pts);

center_=repmat(center,length(pts),1);
dist=sum((pts-center_).^2,2);
maxdist=max(dist);
r=sqrt(maxdist);

end

%% demo
% ptsx=rand(10,1)*10;ptsy=rand(10,1)*10;
% figure;plot(ptsx,ptsy,'b*','MarkerSize',6);axis equal;hold on;
% [center,r]=CircleBox_2D(ptsx,ptsy);
% crv=nrbcirc(r,center);
% nrbctrlplot(crv);
