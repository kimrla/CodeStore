function [tsrf, pta]=tnrbtrimretri(tsrf, p2t2, ed2, pts2, pnts2)

% tnrbtrimretri: Re-triangulation of nodes close to the trimming curve
% 
% Calling Sequences:
% 
%       [tsrf, pta]=tnrbtrimretri(tsrf, p2t2, ed2, pts2, pnts2)
% 
% INPUTS:
% 
%       tsrf - Triangular representation of a nurbs surface.
% 
%       p2t2 - The relations from points to triangles of tri-nurbs 
%                  surface 1. See also tnrbpts2tri.
%
%       ed2  -  Edges of the triangles of the tri-nurbs surface that
%                   intersected with the curve.
%
%       pts2 - Parametric points of the surface.
%
%       pnts2 - The nearest points of the surface to given points.
%
% OUTPUT: 
%
%       tsrf - The tri-nurbs surface after re-triangulation of nodes 
%                close to the trimming curve.
%
%       pta - Indexes of the poits added to the tri-nurbs surface.
%

% Get a band on the tri-nubrs surface for the trimming curve 
tri2=tnrbtrimband(tsrf, p2t2, ed2);

% Get the boundary of the trimming band of a tri-nurbs surface
bp2=tnrbboundary(tsrf, p2t2, tri2);

% Remove triangles and nodes for a tri-nurbs surface
[tsrf, bp2]=tnrbdeletetri(tsrf, tri2, bp2);

% Triangulation of nodes on the trimming band
nt=0; 
nd=zeros(length(bp2),1);
for i=1:length(bp2)
    nd(i)=length(bp2{i});
    nt=nt+nd(i)-1;
end
nodes2=zeros(nt,2);
points2=zeros(nt,3);
pt2=zeros(nt,1);
t=0; 
for i=1:length(bp2)
    nodes2(t+1:t+nd(i)-1,:)=tsrf.nodes(bp2{i}(1:nd(i)-1),:);
    points2(t+1:t+nd(i)-1,:)=tsrf.points(bp2{i}(1:nd(i)-1),:);
    pt2(t+1:t+nd(i)-1)=bp2{i}(1:nd(i)-1);
    bp2{i}=[(t+1:t+nd(i)-1)'; t+1];
    t=t+nd(i)-1;
end
dc=norm(pts2(1,:)-pts2(end,:));
np=length(pts2);
if dc<tsrf.seeds(2)*1e-6;    
    nodes1=pts2(1:np-1,:);
    points1=pnts2(1:np-1,:);
    pt1=[nt+(1:np-1)'; nt+1];
else
    nodes1=pts2(1:np,:);
    points1=pnts2(1:np,:);
    pt1=nt+(1:np)';
end
nodes=[nodes2; nodes1];
points=[points2; points1];
trib=delaunay(nodes(:,1), nodes(:,2));

% Check the direction of the curves and direct them
if length(bp2)==2
    pt1=crvdirect(nodes, pt1);
    t1=nodes(pt1(1),:);
    [~, t21]=nearpnts(t1, nodes(bp2{1},:));
    t21=nodes(bp2{1}(t21),:);
    [~, t22]=nearpnts(t1, nodes(bp2{2},:));
    t22=nodes(bp2{2}(t22),:);
    dr0=vecnorm(nodes(pt1(2),:)-nodes(pt1(1),:));
    dr1=vecnorm(t21-t1);
    dr2=vecnorm(t22-t1);
    dn=cross([dr1, 0], [dr0, 0]);
    dr=ones(1,2);
    if dn(3)<0
        dr(1)=-1;
    end
    dn=cross([dr2, 0], [dr0, 0]);
    if dn(3)<0
        dr(2)=-1;
    end
    eb=cell(1,length(bp2));
    for i=1:length(bp2)
        bp2{i}=crvdirect(nodes, bp2{i}, dr(i));        
        eb{i}=[bp2{i}(1:end-1), bp2{i}(2:end)];
    end
elseif length(bp2)==1
    bp2{1}=crvdirect(nodes, bp2{1});
    eb{1}=[bp2{1}(1:end-1), bp2{1}(2:end)];
end

% Remove nodes outside the boundary of trimming band
p2t=cell(1, length(nodes));
for i=1:length(trib)
    for j=1:3
        p2t{trib(i,j)}=[p2t{trib(i,j)}, i];
    end
end
rt=true(length(trib),1);
for r=1:length(eb)
    ebr=eb{r};
    sebr=sort(ebr, 2);
    for i=1:length(ebr)
        ei=ebr(i,:);
        tri=tnrbedge2tri(p2t, ei);
        if length(tri)==2
            for j=1:length(tri)
                if rt(tri(j))
                    pj=trib(tri(j), (trib(tri(j),:)~=ei(1)) & (trib(tri(j),:)~=ei(2))); 
                    for k=1:2
                        ek=[pj, ei(k)];
                        dm=min(DistanceMatrix(sort(ek), sebr));
                        if dm~=0
                            dr1=vecnorm(nodes(ei(2),:)-nodes(ei(1),:));
                            pm=(nodes(pj,:)+nodes(ei(k),:))/2;
                            dr2=vecnorm(pm-nodes(ei(1),:));
                            dn=cross([dr2,0], [dr1, 0]);
                            if dn(3)>1e-3
                                rt(tri(j))=false;
                                trik=tnrbedge2tri(p2t, ek);
                                trik=trik(trik~=tri(j));
                                if ~isempty(trik)
                                    dm1=min(DistanceMatrix(sort(trib(trik(1),[1,2])), sebr));
                                    dm2=min(DistanceMatrix(sort(trib(trik(1),[2,3])), sebr));
                                    dm3=min(DistanceMatrix(sort(trib(trik(1),[3,1])), sebr));
                                    if dm1~=0 && dm2~=0 && dm3~=0
                                        rt(trik)=false;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
trib=trib(rt,:);

% Add new triangles to a tri-nurbs surface
ida=RemDuplicate(pt1);
pe=pt2; pta=tsrf.numbers(1)+(1:length(ida))';
pt=[pe; pta]; 
tsrf.numbers(1)=tsrf.numbers(1)+length(ida);
tsrf.numbers(2)=tsrf.numbers(2)+length(trib);
tsrf.nodes=[tsrf.nodes; nodes(ida,:)];
tsrf.points=[tsrf.points; points(ida,:)];
for i=1:length(trib)
    for j=1:3
       trib(i,j)=pt(trib(i,j)); 
    end
end
tsrf.delaunay=[tsrf.delaunay; trib];

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
% % Plot results
% figure; hold on; 
% tnrbplot(tsrf);
% tnrbplot(tcrv);
% axis equal;
% title('The original surface and the trimming curve.');
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
% % Plot results
% figure; hold on; 
% tnrbplot(tsrf);
% tnrbplot(tcrv);
% plot3(tsrf.points(pta,1), tsrf.points(pta,2), tsrf.points(pta,3), 'ro')
% axis equal;
% title('The surface after re-triangulation.');
% 
% figure; hold on;
% triplot(tsrf.delaunay, tsrf.nodes(:,1), tsrf.nodes(:,2)); 
% plot(tsrf.nodes(pta,1), tsrf.nodes(pta,2), 'ro')
% plot(pts2(:,1), pts2(:,2));
% axis equal;
% title('Parametric domian of the surface after re-triangulation.');






