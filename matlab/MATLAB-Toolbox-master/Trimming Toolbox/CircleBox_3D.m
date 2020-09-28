function [ center,r ] = CircleBox_3D( varargin )
% Construct 3D sphere box from scattered points or the coordinates.

% Input:
%     ptx,pty: The x and y coordinates of points, or the scattered points pts.
% Output:
%     center: The centroid point of the sphere box.
%     r: The radius of the circle.


if nargin==3
    pts=[varargin{1}(:),varargin{2}(:),varargin{3}(:)];
else
    pts=varargin{1};
end

center=sum(pts)/length(pts);

center_=repmat(center,length(pts),1);
dist=sum((pts-center_).^2,2);
maxdist=max(dist);
r=sqrt(maxdist);

end

