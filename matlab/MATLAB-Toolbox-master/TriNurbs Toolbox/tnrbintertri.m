function [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri(tnrb1, tnrb2, p2t1, p2t2, p1, p2)

% tnrbintertri: Get the intersection points of two tri-nurbs surfaces
% 
% Calling Sequences:
% 
%       [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri(tnrb1, tnrb2, p2t1, p2t2, p1, p2)
% 
% INPUTS:
% 
%       tnrb1, tnrb2 - Triangular representation of nurbs surface.
% 
%       p2t1, p2t2 - The relations from points to triangles of two 
%                     tri-nurbs surface. See also tnrbpts2tri.
% 
%       p1, p2 - Indexes of nearest points between two surfaces.
%                    See also nearpnts.
%
% OUTPUT: 
%
%      ed1  -  Edges of the triangles of tri-nurbs surface 1 that
%                   intersecte with tri-nurbs surface 2.
%
%      tri2   -  Triangles tri-nurbs surface 2 that intersecte with 
%                  tri-nurbs surface 1.
%
%       pts1, pts2 - Parametric intersection points of the two surfaces.
%
%       pnts1, pnts2 - Intersection points of the two surfaces.
%

% Get approximated intersections by triangles
nd=length(p1); t=0;
for k=1:nd 
    edges=tnrbpt2edges(tnrb1, p2t1, p1(k));
    m=length(edges);
    n=length(p2t2{p2(k)});
    for i=1:m
        for j=1:n
            ed=[p1(k), edges(i)];
            tri=tnrb2.delaunay(p2t2{p2(k)}(j),:);
            [tpnts2, tpnts1, x]=tnrbinterline2tri(tnrb1, tnrb2, ed, tri);
            tr=tnrb2.nodes(tri, :);
            line=tnrb1.nodes(ed, :); 
            if ~isempty(x)
                tpts1x=(1-x(3,:))*line(1,1)+x(3,:)*line(2,1);%x是线段的局部坐标，tpts1是线段落在曲面参数域内的整体坐标
                tpts1y=(1-x(3,:))*line(1,2)+x(3,:)*line(2,2);
                tpts1=[tpts1x(:), tpts1y(:)];
                tpts2=tripoint(tr(1,:), tr(2,:), tr(3,:), x(1,:), x(2,:));
                tk=size(x, 2);
                ed1(t+1:t+tk,:)=ed;
                tri2(t+1:t+tk,:)=tri;
                pts1(t+1:t+tk,:)=tpts1;
                pts2(t+1:t+tk,:)=tpts2;
                pnts1(t+1:t+tk,:)=tpnts1;
                pnts2(t+1:t+tk,:)=tpnts2;
                t=t+tk;
            end
        end
    end
end
ed1=sort(ed1, 2);
[ed1, id]=RemDuplicate(ed1);
tri2=tri2(id,:);
pts1=pts1(id,:);
pts2=pts2(id,:);
pnts1=pnts1(id,:);
pnts2=pnts2(id,:);


%% demo
% % The mesh seed length (h0)
% h0=0.9;
% 
% % Create a nurbs sphere
% center=[5,1,4];
% circ=nrbcirc(4, center, 0, pi);
% srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
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
% % Get the relations from points to triangles of tri-nurbs
% p2t1=tnrbpts2tri(tnrb1);
% p2t2=tnrbpts2tri(tnrb2);
% 
% % Intersections of two tri-nurbs surfaces
% [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri(tnrb1, tnrb2, p2t1, p2t2, p1, p2);
% 
% % Plot results
% figure; hold on; 
% tnrbplot(tnrb1); 
% tnrbplot(tnrb2);
% plot3(pnts1(:,1), pnts1(:,2), pnts1(:,3), 'ro'); 
% plot3(pnts2(:,1), pnts2(:,2), pnts2(:,3), 'r*'); 
% axis equal; view(3);
% title('Geometric grid');
% 
% figure; hold on;
% triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2)); 
% plot(pts1(:,1), pts1(:,2), 'k.', 'MarkerSize', 13); 
% title('Parametric mesh of surface 1'); 
% axis equal;
% 
% figure; hold on;
% triplot(tnrb2.delaunay, tnrb2.nodes(:,1), tnrb2.nodes(:,2)); 
% plot(pts2(:,1), pts2(:,2), 'k.', 'MarkerSize', 13); 
% title('Parametric mesh of surface 2'); 
% axis equal;


