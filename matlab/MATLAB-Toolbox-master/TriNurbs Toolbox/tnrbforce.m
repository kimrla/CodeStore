function F=tnrbforce(tsrf, p2t2, pt, h0)

% tnrbforce: Get the force at a point on a tri-nurbs surfaces
% 
% Calling Sequences:
% 
%       F=tnrbforce(tsrf, p2t2, pt, h0)
% 
% INPUTS:
% 
%       tsrf - Triangular representation of a nurbs surface.
% 
%       p2t2 - The relations from points to triangles of tri-nurbs 
%                  surface 1. See also tnrbpts2tri.
% 
%       pt - Index of a point of the tri-nurbs surface.
%
%       h0 - Mesh seed length.
% 
% OUTPUT: 
% 
%       F  -  Force at the point on the tri-nurbs surfaces computed 
%              from length of the edges passed through the point.
%

F=zeros(1,3);
eds=tnrbpt2edges(tsrf, p2t2, pt);
for j=1:length(eds)
    dr=tsrf.points(pt,:) - tsrf.points(eds(j),:);
    hj=norm(dr);
    f=h0-hj;
    F=F+f*dr/hj;
end

%% demo
% % The mesh seed length (h0) and force facter (k0)
% h0=0.5; k0=1;
% 
% % Create a plane square and a plane cuve
% crv=nrbcirc(2, [5,5], 0, 2*pi);
% lin1=nrbline([0,0], [10,0]);
% lin2=nrbline([0,10], [10,10]);
% srf=nrbruled(lin1, lin2);
% 
% % Transform a nurbs surface into triangular representation
% tcrv=nrb2tri(crv, 0.8*h0);
% tsrf=nrb2tri(srf, h0);
% 
% % The nearest points from the surface to the curve
% tol=max([tsrf.seeds(1), tcrv.seeds(1)]);
% [p1, p2, d]=nearpnts(tsrf.points, tcrv.points, tol);
% 
% % Get the relations from points to triangles of tri-nurbs
% p2t2=tnrbpts2tri(tsrf);
% 
% % Get the intersection points of a curve with a tri-nurbs surface
% pnts1=tcrv.points;
% [pts2, pnts2, dist]=tnrbpoint2srf(tsrf, tcrv.points);
% ed2=tnrbinterline(tsrf, tcrv, p2t2, p1, p2);
% 
% % Re-triangulation of nodes close to the trimming curve
% [tsrf, pta]=tnrbtrimretri(tsrf, p2t2, ed2, pts2, pnts2);
% 
% % Get the boundary of the tri-nurbs surface
% p2t2=tnrbpts2tri(tsrf);
% bp=triboundary(tsrf, p2t2);
% bp=bp{1};
% 
% % Optimize the triangles
% pfix=sort(RemDuplicate([pta; bp]));
% pop=true(tsrf.numbers(1),1);
% pop(pfix)=false;
% pop=sort(find(pop));
% np=length(pop);
% F1=zeros(np, 3);
% V=F1;
% for i=1:np
%     F1(i,:)=tnrbforce(tsrf, p2t2, pop(i), h0);
%     
% end
% 
% % Plot results
% figure; hold on; 
% tnrbplot(tsrf);
% tnrbplot(tcrv);
% plot3(tsrf.points(pta,1), tsrf.points(pta,2), tsrf.points(pta,3), 'ro')
% plot3(tsrf.points(bp,1), tsrf.points(bp,2), tsrf.points(bp,3), 'r')
% quiver(tsrf.points(pop,1), tsrf.points(pop,2), F1(:,1), F1(:,2));
% axis equal;
% title('The surface after re-triangulation.');



