function [tsrf, bp]=tnrbdeletetri(tsrf, tri, bp)

% tnrbdeletetri: Remove triangles and nodes for a tri-nurbs surface
% 
% Calling Sequences:
% 
%       tsrf=tnrbdeletetri(tsrf, tri)
% 
%       [tsrf, bp]=tnrbdeletetri(tsrf, tri, bp)
% 
% INPUTS:
% 
%       tsrf - A tri-nurbs surface.
% 
%       tri - Indexes of triangles to be removed.
%
%       bp  -  Indexes of nodes need to be renumber after the deletion.
%
% OUTPUT: 
%
%      tsrf - The tri-nurbs surface after deletion.
%
%      bp - Indexes of nodes after renumberinmg.
%

% Delete triangles
rv=true(tsrf.numbers(2),1); 
rv(tri)=false;
tsrf.numbers(2)=tsrf.numbers(2)-length(tri);
tsrf.delaunay=tsrf.delaunay(rv,:);

% Delete and renumber points and nodes
nd=tsrf.numbers(1);
vd=false(nd,1);
vv=zeros(nd,1);
for i=1:tsrf.numbers(2)
    for j=1:3
        vd(tsrf.delaunay(i,j))=true;
    end
end
vm=find(vd);
tsrf.numbers(1)=length(vm);
vv(vm)=1:tsrf.numbers(1);
for i=1:tsrf.numbers(2)
    for j=1:3
        tsrf.delaunay(i,j)=vv(tsrf.delaunay(i,j));
    end
end
if nargin==3
    for i=1:length(bp)
        for j=1:length(bp{i})
            bp{i}(j)=vv(bp{i}(j));
        end
    end
end
tsrf.nodes=tsrf.nodes(vd,:);
tsrf.points=tsrf.points(vd,:);

if nargin==2
    bp=[];
end

%% demo
% % The mesh seed length (h0)
% h0=0.5;
% 
% % Create a plane square and a plane cuve
% crv=nrbcirc(2, [5,5], 0, 2*pi);
% lin1=nrbline([0,0], [10,0]);
% lin2=nrbline([0,10], [10,10]);
% srf=nrbruled(lin1, lin2);
% 
% % Transform a nurbs surface into triangular representation
% tcrv=nrb2tri(crv, h0);
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
% % Get a band on the tr-nubrs surface for the trimming curve 
% tri2=tnrbtrimband(tsrf, p2t2, ed2);
% 
% % Get the boundary of trimming band of a tri-nurbs surface
% bp2=tnrbboundary(tsrf, p2t2, tri2);
% 
% % Remove triangles and nodes for a tri-nurbs surface
% [tsrf, bp2]=tnrbdeletetri(tsrf, tri2, bp2);
% 
% The nodes to be retained
% fixnd=[]; 
% for i=1:length(bp2)
%     fixnd=[fixnd; bp2{i}];
% end
% fixnd=RemDuplicate(fixnd);

% % Plot results
% figure; hold on; 
% tnrbplot(tsrf);
% tnrbplot(tcrv);
% plot(pnts1(:,1), pnts1(:,2), 'ro');
% plot(pnts2(:,1), pnts2(:,2), 'r*');
% axis equal;

% figure; hold on;
% triplot(tsrf.delaunay, tsrf.nodes(:,1), tsrf.nodes(:,2)); 
% plot(pts2(:,1), pts2(:,2), 'k.', 'MarkerSize', 13); 
% for i=1:length(bp2)
%     plot(tsrf.nodes(bp2{i},1), tsrf.nodes(bp2{i},2), 'k', 'LineWidth', 1); 
% end
% axis equal;


