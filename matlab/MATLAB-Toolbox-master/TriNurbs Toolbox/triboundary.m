function bp=triboundary(tsrf, p2t)

% triboundary: Get the boundary of a tnrb surface.
% 和tnrbboundary比较（获得给定三角形的边界，和tnrbtrimband配合使用）
% Calling Sequences:
% 
%     bp=tnrbboundary(tsrf, p2t)
% 
% INPUTS:
%
%      tsrf - Triangular representation of a nurbs (tri-nurbs) surface.
%
%      p2t - The relations from points to triangles. See tnrbpts2tri.
%
% OUTPUT:
% 
%    bp - A cell array. Each cell contains a vector of boundary nodes.
%  

% Find edges with only one triangle
tris=(1:tsrf.numbers(2))';
e1=tsrf.delaunay(tris, [1,2]);
e2=tsrf.delaunay(tris, [2,3]);
e3=tsrf.delaunay(tris, [3,1]);
et=sort([e1; e2; e3], 2);
et=RemDuplicate(et);
eb=[];
for i=1:length(et)
    te1=p2t{et(i,1)};
    te2=p2t{et(i,2)};
    tri=[te1, te2]';
    [~, id]=RemDuplicate(tri);
    tri=tri(~id);
    if length(tri)==1
        eb=[eb; i];
    end
end
et=et(eb,:);

% Sort the boundary edges
r=1; t=2;
k(1:2,1)=et(1,:);
pp=[2, 1];
ne=length(et);
id=true(ne,1);
id(1)=false;
di=zeros(ne,1);
di(1)=1;
for i=2:ne
    p1=find(et(:,1)==k(t));
    p2=find(et(:,2)==k(t));
    j1=1*ones(size(p1));
    j2=2*ones(size(p2));
    p=[p1;p2];
    j=[j1; j2];
    q=find(p~=di(i-1));
    if ~isempty(q) && id(p(q(1)))
        k(t+1)=et(p(q(1)), pp(j(q(1))));
        id(p(q(1)))=false;
        t=t+1;
        di(i)=p(q(1));
    end
    if k(t)==k(1)
        bp{r}=k;
        r=r+1;
        k=find(id);
        if isempty(k)
            break;
        end
        k(1:2,1)=et(k(1),:);
        t=2;
    end
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
% % Plot results
% figure; hold on; 
% tnrbplot(tsrf);
% tnrbplot(tcrv);
% plot3(tsrf.points(pta,1), tsrf.points(pta,2), tsrf.points(pta,3), 'ro')
% plot3(tsrf.points(bp,1), tsrf.points(bp,2), tsrf.points(bp,3), 'r')
% axis equal;
% title('The surface after re-triangulation.');
% 
% figure; hold on;
% triplot(tsrf.delaunay, tsrf.nodes(:,1), tsrf.nodes(:,2)); 
% plot(tsrf.nodes(pta,1), tsrf.nodes(pta,2), 'ro')
% plot(pts2(:,1), pts2(:,2));
% axis equal;
% title('Parametric domian of the surface after re-triangulation.');







