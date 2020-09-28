function tr=uniqpt2tri(qnrb, pk)

% uniqpt2tri: Get all triangles related with a unique point on a quad-nurbs.
% 
% Calling Sequences:
% 
%     tri=uniqpt2tri(qnrb, pk)
% 
% INPUTS:
%
%      qnrb - Quadrangular representation of a nurbs (quad-nurbs) surface.
%
%      pk  -  A index of a unqiue point on a quad-nurbs surface.
%
% OUTPUT:
%
%    tri - The triangles related with the point.
%  

if pk>size(qnrb.UniqPt, 1)
    error('Index of the point is larger than the total number of unique points!');
end

dump=qnrb.UniqId{2};
pt2tri=qnrb.pt2tri;
np=length(dump{pk});
tr=[];
for i=1:np
    tr=[tr, pt2tri{dump{pk}(i)}];
end
tr=unique(tr);


%% demo
% % The mesh seed length (h0)
% h0=1.5;
% 
% % Create a nurbs sphere
% L=[7,10,0; 5,1,10];
% center=[5,5,4];
% circ=nrbcirc(4, center, 0, pi);
% srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% srf=srf1;
% 
% % Transform a nurbs surface into quadrangular representation
% qnrb=nrb2quad(srf, h0);
% 
% % Get unique edges of quadrangles related to a unique point
% k=120;
% tr=uniqpt2tri(qnrb, k);
% 
% % Plot the surface and related quadrangle and line
% figure; hold on;
% trisurf(qnrb.quad, qnrb.points(:,1), qnrb.points(:,2), qnrb.points(:,3));
% trimesh(qnrb.tri(tr,:), qnrb.points(:,1), qnrb.points(:,2), qnrb.points(:,3));
% view(3); axis equal;







