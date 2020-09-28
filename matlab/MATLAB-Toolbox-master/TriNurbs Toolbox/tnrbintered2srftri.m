function [de, p1, pt1, p2, pt2, et]=tnrbintered2srftri(tnrb1, tnrb2, ed1, pts2, sqt, tri, it, n)

% tnrbintered2srftri: Get the intersection points of two tri-nurbs surfaces through a known edge 
%                           and a given triangle that includes the edge
% 
% Calling Sequences:
% 
%       [de, p1, pt1, p2, pt2, et]=tnrbintered2srftri(tnrb1, tnrb2, ed1, pts2, sqt, tri)
% 
%       [de, p1, pt1, p2, pt2, et]=tnrbintered2srftri(tnrb1, tnrb2, ed1, pts2, sqt, tri, it)
% 
%       [de, p1, pt1, p2, pt2, et]=tnrbintered2srftri(tnrb1, tnrb2, ed1, pts2, sqt, tri, [], n)
% 
%       [de, p1, pt1, p2, pt2, et]=tnrbintered2srftri(tnrb1, tnrb2, ed1, pts2, sqt, tri, it, n)
% 
% INPUTS:
% 
%       tnrb1, tnrb2 - Triangular representation of two nurbs surface.
%
%      ed1  -  Edges of the triangles of tri-nurbs surface 1 that
%                   intersected with tri-nurbs surface 2.
% 
%       pts2 - Approximated parametric intersections of the two surfaces.
%
%       sqt  -  An index of ed1. 
%
%       tri  -   The triangle that includes the edge ed1(sqt)
%
%       it -  The number of iterations.
%
%       n  -  The number of grids in each subdomain used for surface.
%
% OUTPUT: 
%
%       de  -  Distances between the edges of the triangle of tnrb1 and the surface tnrb2.
%
%       p1, p2 - Parametric intersection points of the two surfaces.
%
%       pt1, pt2 - Intersection points of the two surfaces.
%
% Discription:
%
%       See also tnrbintersort.
%

if nargin==6
    it=10; n=8; 
elseif nargin==7
    n=8;
end
if isempty(it)
    it=10;
end

td=tnrb1.delaunay(tri,:);
et=ed1(sqt,:);
pe=td(((td~=et(1)) & (td~=et(2))));
if length(pe)>1
    error('The triangle does not include the edge!');
end
et=sort([pe, et(1); pe, et(2)], 2);
p1=(tnrb1.nodes(et(:,1),:)+tnrb1.nodes(et(:,2),:))/2;
p2=[pts2(sqt, :); pts2(sqt, :)];
pt1=zeros(2,3); 
pt2=pt1;
de=zeros(2,1);
for i=1:2
    [de(i), p1(i,:), pt1(i,:), p2(i,:), pt2(i,:)]=tnrbintered2srf(tnrb1, tnrb2, et(i,:), p1(i,:), p2(i,:), it, n);
end




