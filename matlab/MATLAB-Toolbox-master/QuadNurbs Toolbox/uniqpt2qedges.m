function [edges, qd]=uniqpt2qedges(qnrb, pk)

% uniqpt2qedges: Get all edges related with a unique point on a quad-nurbs.
% 
% Calling Sequences:
% 
%     [edges, qd]=uniqpt2qedges(qnrb, pk)
% 
% INPUTS:
%
%      qnrb - Quadrangular representation of a nurbs (quad-nurbs) surface.
%
%      pk  -  Indexes of unqiue points on a quad-nurbs surface.
%
% OUTPUT:
% 
%    edges - Point indexes of all edges related with the points 
%                 on the quad-nurbs surface.
%
%    qd - The quadrangles related with the point.
%  

if max(pk)>size(qnrb.UniqPt, 1)
    error('Index of the point is larger than the total number of unique points!');
end

points=qnrb.points;
dump=qnrb.UniqId{2};
pt2qd=qnrb.pt2qd;
h0=qnrb.seeds;
nk=length(pk);
qd=zeros(1, 4*nk);
t=0;
for j=1:nk
    np=length(dump{pk(j)});
    for i=1:np
        tk=length(pt2qd{dump{pk(j)}(i)});
        qd(t+1:t+tk)=pt2qd{dump{pk(j)}(i)};
        t=t+tk;
    end
end
qd(t+1:end)=[];
qd=unique(qd);
edges=[qnrb.quad(qd, [1,2])
    qnrb.quad(qd, [2,3])
    qnrb.quad(qd, [3,4])
    qnrb.quad(qd, [1,4])];
edges=RemDuplicate(sort(edges, 2));
d=points(edges(:,1),:)-points(edges(:,2),:);
d=sqrt(sum(d.^2, 2));
p=d>h0*1e-6;
edges=edges(p,:);


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
% k=26;
% [edges, qd]=uniqpt2qedges(qnrb, k);
% 
% % Plot the surface and related quadrangle and line
% figure; hold on;
% trisurf(qnrb.quad, qnrb.points(:,1), qnrb.points(:,2), qnrb.points(:,3));
% trimesh(qnrb.quad(qd,:), qnrb.points(:,1), qnrb.points(:,2), qnrb.points(:,3));
% for i=1:size(edges, 1)
%     plot3(qnrb.points(edges(i,:),1), qnrb.points(edges(i,:),2), qnrb.points(edges(i,:),3), 'r');
% end
% view(3); axis equal;







