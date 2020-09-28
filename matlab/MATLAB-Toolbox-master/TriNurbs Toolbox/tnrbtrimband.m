function [tri, pnts]=tnrbtrimband(tsrf, p2t, ed)

% tnrbtrimband: Get a band on the tr-nubrs surface for the trimming curve 
% 
% Calling Sequences:
% 
%       tri=tnrbtrimband(tsrf, p2t, ed)
% 
%       [tri, pnts]=tnrbtrimband(tsrf, p2t, ed)
% 
% INPUTS:
% 
%       tsrf - A tri-nurbs surface.
% 
%       p2t - The relations from points to triangles of the surface.
%
%      ed1  -  Edges of the triangles of the tri-nurbs surface that
%                   intersected with a curve or surface.
%
% OUTPUT: 
%
%      tri  -  Tiangles of the band.
%
%      pnts - Vertexes of the triangles.
%

tri=[];
for i=1:length(ed)
    tri=[tri, p2t{ed(i,1)}, p2t{ed(i,2)}];
end
tri=RemDuplicate(tri');

if nargout==2
    pnts=[];
    for i=1:length(tri)
        pnts=[pnts, tsrf.delaunay(tri(i),:)];
    end
    pnts=RemDuplicate(pnts');
end

