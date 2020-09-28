function p2t=tnrbpts2tri(tnrb)

% tnrbpts2tri: Get the relations from points to triangles of tri-nurbs.
% 
% Calling Sequences:
% 
%     p2t=tnrbpts2tri(tnrb)
% 
% INPUTS:
%
%      tnrb - Triangular representation of a nurbs surface.
%
% OUTPUT:
% 
%    p2t - The relations from points to triangles.
%  

nd=tnrb.numbers(1);
nt=tnrb.numbers(2);
tri=tnrb.delaunay;
p2t=cell(1, nd);
for i=1:nt
    for j=1:3
        p2t{tri(i,j)}=[p2t{tri(i,j)}, i];
    end
end


%% demo
% % The mesh seed length (h0)
% h0=0.2;
% 
% % Create a nurbs sphere
% circ=nrbcirc(1, [0,0], 0, pi);
% srf=nrbrevolve(circ, [0,0,0], [1,0,0], 2*pi);
% 
% % Transform a nurbs surface into triangular representation
% tnrb=nrb2tri(srf, h0);
% 
% % Plot the triangular nurbs surface
% figure; hold on;
% trisurf(tnrb.delaunay, tnrb.points(:,1), tnrb.points(:,2), tnrb.points(:,3));
% axis equal; view(3);
% 
% % The relations from points to triangles
% p2t=tnrbpts2tri(tnrb);
% 
% % Plot a node and related triangles
% k=219;
% tri=tnrb.delaunay;
% trik=tri(p2t{k}, :);
% trimesh(trik, tnrb.points(:,1), tnrb.points(:,2), tnrb.points(:,3));
% plot3(tnrb.points(k,1), tnrb.points(k,2), tnrb.points(k,3), 'ro')






