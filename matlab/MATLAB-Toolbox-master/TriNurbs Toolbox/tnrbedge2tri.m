function tri=tnrbedge2tri(p2t, ed)

% tnrbedge2tri: Get the relations from an edge to triangles.
% 
% Calling Sequences:
% 
%     tri=tnrbedge2tri(tnrb, ed)
% 
% INPUTS:
%
%      p2t - The relations from points to triangles.
%
%      ed  -  An edge.
%
% OUTPUT:
% 
%    tri - The triangle has the edge.
%  

tri1=p2t{ed(1)};
tri2=p2t{ed(2)};
tri=[tri1, tri2];
[~, id]=RemDuplicate(tri');
tri=tri(~id);

%% demo
% % The mesh seed length (h0)
% h0=1.5;
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
% % Get the intersection points of two tri-nurbs surfaces
% [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri(tnrb1, tnrb2, p2t1, p2t2, p1, p2);
% 
% % The relations of an edge with a triangle
% k=1; ed=ed1(k,:);
% tri=tnrbedge2tri(p2t1, ed1(k,:)); 
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
% plot(tnrb1.nodes(ei,1), tnrb1.nodes(ei,2), 'g.', 'MarkerSize', 15);
% plot(tnrb1.nodes(ei,1), tnrb1.nodes(ei,2), 'g', 'LineWidth', 1.5);
% plot(pts1(:,1), pts1(:,2), 'r*');
% title('Parametric mesh of surface 1'); 
% axis equal;
% 
% figure; hold on;
% triplot(tnrb2.delaunay, tnrb2.nodes(:,1), tnrb2.nodes(:,2)); 
% plot(pts2(:,1), pts2(:,2), 'r*');
% title('Parametric mesh of surface 2'); 
% axis equal;



