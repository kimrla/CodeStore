function [m,linid]=findtri(tri,pnt,st1,st2)
% Find the triangulars the line containing one intersection point belongs
%       to. If the intersection is in the four edge of parameter
%       domain,then return 1 triangular ID. If the intersection is
%       generally in other grid lines, then return 2 triangular IDs.
% Input:
%   tri: triangulation structure, see Triangulation or delaunayTriangulation
%   pnt: coordinates of the point in parameter domain
%   st1: knots vector in x direction
%   st2: knots vector in y direction
% output:
%   m: ID of the triangular, or 2 IDs of the 2 triangulars, which
%       corresponding to the tri.ConnectivityList
%   linid: ID of the vertices in one line, which is a 2-element vector

TR=tri.ConnectivityList;
pnts=tri.Points;
TR=sort(TR,2);

linid=findvert(pnt,pnts,st1,st2);
linid=sort(linid);

% if (pnts(linid(1),1)==pnts(linid(2),1))
% if (nargout==1)
m=find( (TR(:,1)==linid(1) & TR(:,2)==linid(2)) | (TR(:,1)==linid(1) & TR(:,3)==linid(2)) | (TR(:,2)==linid(1) & TR(:,3)==linid(2)) ); 
end

    


    
    

% 同样不考虑交点恰好位于网格线交点的情况





