function edges=tnrbpt2edges(tnrb, p2t, pt)

% tnrbpt2edges: Get all edges related with a point on a tri-nurbs.
% 
% Calling Sequences:
% 
%     edges=tnrbpt2edges(tnrb, p2t, pt)
% 
% INPUTS:
%
%      tnrb - Triangular representation of a nurbs (tri-nurbs) surface.
%
%      p2t - The relations from points to triangles. See tnrbpts2tri.
%
%      pt  -  A index of a point on a tri-nurbs surface.
%
% OUTPUT:
% 
%    edges - Point indexes of all edges related with the point 
%                 on the tri-nurbs surface.
%  

tri=p2t{pt};
n=length(tri);
edges=zeros(2*n, 1);
for i=1:n
    p=tnrb.delaunay(tri(i),:)~=pt;
    edges(2*i-1:2*i)=tnrb.delaunay(tri(i),p);
end
edges=RemDuplicate(edges);

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
% % Get the relations from points to edges of tri-nurbs.
% k=58; pt=p1(k);
% p2t=tnrbpts2tri(tnrb1);
% edges=tnrbpt2edges(tnrb1, p2t, pt);
% 
% % Plot the results
% figure; hold on;
% tnrbplot(tnrb1);
% tnrbplot(tnrb2);
% axis equal; view(3);
% title('Geometric grid');
% plot3(tnrb1.points(p1,1), tnrb1.points(p1,2), tnrb1.points(p1,3), 'ro');
% plot3(tnrb2.points(p2,1), tnrb2.points(p2,2), tnrb2.points(p2,3), 'r*');
% 
% for i=1:length(edges)
%     ei=[pt, edges(i)];
%     plot3(tnrb1.points(ei,1), tnrb1.points(ei,2), tnrb1.points(ei,3), 'g.', 'MarkerSize', 15);
%     plot3(tnrb1.points(ei,1), tnrb1.points(ei,2), tnrb1.points(ei,3), 'g', 'LineWidth', 1.5);
% end






