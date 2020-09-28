function t2e=tnrbtri2edge(tnrb, ed, tri)

% tnrbtri2dges: Get the relations from triangles to edges.
% 
% Calling Sequences:
% 
%     t2e=tnrbtri3dges(tnrb, ed, tri)
% 
% INPUTS:
%
%      tnrb - Triangular representation of a nurbs surface.
%
%      ed  -  Edges.
%
%      tri - Triangles inlude edges in ed.
%
% OUTPUT:
% 
%      t2e - The edges related with the triangle.
%  

if isempty(tri)
    error('No triangle is given.');
end
e1=tnrb.delaunay(tri, [1,2]);
e2=tnrb.delaunay(tri, [2,3]);
e3=tnrb.delaunay(tri, [3,1]);
et=sort([e1; e2; e3], 2);
d=DistanceMatrix(et, ed);
[dm, id]=min(d, [], 2);
t2e=id(dm==0);
t2e=RemDuplicate(t2e);

%% demo
% % The mesh seed length (h0)
% h0=0.8;
% 
% % Create a nurbs sphere
% circ=nrbcirc(4, [5,5,4], 0, pi);
% srf1=nrbrevolve(circ, [5,5,4], [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% 
% % Transform a nurbs surface into triangular representation
% tnrb1=nrb2tri(srf1, h0);
% tnrb2=nrb2tri(srf2, h0);
% 
% % The nearest points of the two surfaces
% tol=max([tnrb1.seeds(1), tnrb2.seeds(1)]);
% [p1, p2, d]=nearpnts(tnrb1.points, tnrb2.points, tol);
% 
% % Get the edges of a tri-nurbs surface that intersected with another tri-nurbs surface
% p2t1=tnrbpts2tri(tnrb1);
% p2t2=tnrbpts2tri(tnrb2);
% 
% % Interections of two tri-nurbs surfaces
% [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri(tnrb1, tnrb2, p2t1, p2t2, p1, p2);
% 
% % The relations of an edge with a triangle
% k=1;
% p2t1=tnrbpts2tri(tnrb1);
% tri=tnrbedge2tri(p2t1, ed1(k,:)); 
% 
% % Get the relations from triangles to edges
% t2e=tnrbtri2dges(tnrb1, ed1, tri);
% 
% % Plot the results
% figure; hold on;
% tnrbplot(tnrb1);
% tnrbplot(tnrb2);
% axis equal; view(3);
% title('Geometric grid');
% plot3(pnts1(:,1), pnts1(:,2), pnts1(:,3), 'ro');
% plot3(pnts2(:,1), pnts2(:,2), pnts2(:,3), 'r*');
% ei=ed1(k,:);
% plot3(tnrb1.points(ei,1), tnrb1.points(ei,2), tnrb1.points(ei,3), 'g.', 'MarkerSize', 15);
% plot3(tnrb1.points(ei,1), tnrb1.points(ei,2), tnrb1.points(ei,3), 'g', 'LineWidth', 1.5);
% 
% figure; hold on;
% triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2)); 
% trimesh(tnrb1.delaunay(tri,:), tnrb1.nodes(:,1), tnrb1.nodes(:,2), 'r', 'LineWidth',2); 
% for i=1:length(t2e)
%     ei=ed1(t2e(i),:);
%     plot(tnrb1.nodes(ei,1), tnrb1.nodes(ei,2), 'g.', 'MarkerSize', 15);
%     plot(tnrb1.nodes(ei,1), tnrb1.nodes(ei,2), 'g', 'LineWidth', 1.5);
% end
% plot(pts1(:,1), pts1(:,2), 'r*');
% title('Parametric mesh of surface 1'); 
% axis equal;
% 
% figure; hold on;
% triplot(tnrb2.delaunay, tnrb2.nodes(:,1), tnrb2.nodes(:,2)); 
% plot(pts2(:,1), pts2(:,2), 'r*');
% title('Parametric mesh of surface 2'); 
% axis equal;




