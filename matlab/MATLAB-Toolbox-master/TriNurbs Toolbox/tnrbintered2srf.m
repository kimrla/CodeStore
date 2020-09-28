function [d, pts1, pnts1, pts2, pnts2]=tnrbintered2srf(tnrb1, tnrb2, ed, pts1, pts2, it, n)

% Calling Sequences:
% 
%       [d, pts1, pnts1, pts2, pnts2]=tnrbintered2srf(tnrb1, tnrb2, ed, pts1, pts2)
% 
%       [d, pts1, pnts1, pts2, pnts2]=tnrbintered2srf(tnrb1, tnrb2, ed, pts1, pts2, it)
% 
%       [d, pts1, pnts1, pts2, pnts2]=tnrbintered2srf(tnrb1, tnrb2, ed, pts1, pts2, it, n)
% 
% INPUTS:
% 
%       tnrb1, tnrb2 - Triangular representation of two nurbs surface.
%
%       ed -  Indexes of an edge on tnrb1.
% 
%       pts1, pts2 - Approximated parametric intersections of the two surfaces.
%
%       it -  The number of iterations.
%
%       n  -  The number of grids in each subdomain used for surface.
%
% OUTPUT: 
%
%       d  -  Distances between the edge of tnrb1 and the surface tnrb2.
%
%       pts1, pts2 - Parametric intersection points of the two surfaces.
%
%       pnts1, pnts2 - Intersection points of the two surfaces.
%

if nargin==5
    it=8; n=8; 
elseif nargin==6
    n=8;
end
if isempty(it)
    it=8;
end

edgseed = @(tnrb, r, ei, k) (1-r)*tnrb.nodes(ei(1), k)+r*tnrb.nodes(ei(2), k);
u1=tnrb1.nodes(ed(1), 1);
u2=tnrb1.nodes(ed(2), 1);
v1=tnrb1.nodes(ed(1), 2);
v2=tnrb1.nodes(ed(2), 2);
du=u2-u1; dv=v2-v1; 
if du~=0
    r1=(pts1(1)-u1)/du;
else
    r1=(pts1(2)-v1)/dv;
end
pts1=zeros(2, n);
dr=0.5; 
ds=tnrb2.seeds(2); 
for i=1:it
    r=linspace(max([0, r1-dr]), min([1, r1+dr]), n);
    s=linspace(max([0, pts2(1)-ds]), min([1, pts2(1)+ds]), n);
    t=linspace(max([0, pts2(2)-ds]), min([1, pts2(2)+ds]), n);
    pts1(1,:)=edgseed(tnrb1, r, ed, 1);
    pts1(2,:)=edgseed(tnrb1, r, ed, 2);
    pt1=nrbeval(tnrb1.nurbs, pts1);
    pt2=nrbeval(tnrb2.nurbs, {s, t});
    pt1=pt1'; pt2=pt2(:,:)'; 
    [~, pti, di]=nearpnts(pt1, pt2);
    [d, id]=min(di); 
    pnts1=pt1(id, :); 
    pnts2=pt2(pti(id), :); 
    r1=r(id); 
    i2=rem(pti(id), n); 
    j2=fix(pti(id)/n)+1; 
    if i2==0
        j2=j2-1; i2=n;
    end
    pts2=[s(i2), t(j2)];
    dr=dr/2; ds=ds/2;
end
pts1=pts1(:,id)';


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
% % Plot the results
% figure; hold on;
% tnrbplot(tnrb1);
% tnrbplot(tnrb2);
% axis equal; view(3);
% title('Geometric grid');
% plot3(tnrb1.points(p1,1), tnrb1.points(p1,2), tnrb1.points(p1,3), 'ro');
% plot3(tnrb2.points(p2,1), tnrb2.points(p2,2), tnrb2.points(p2,3), 'r*');
% 
% % Get the relations from points to edges of tri-nurbs.
% k=5; i=1; 
% pt1=p1(k); pt2=p2(k); 
% p2t=tnrbpts2tri(tnrb1);
% edges=tnrbpt2edges(tnrb1, p2t, pt1);
% ed=[pt1, edges(i)];
% 
% % Get intersections of a edge of a tri-nurbs surface with another tri-nurbs surface
% pts1=(tnrb1.nodes(ed(1),:)+tnrb1.nodes(ed(2),:))/2;
% pts2=tnrb2.nodes(pt2, :);
% [dt, pts1, pnts1, pts2, pnts2]=tnrbintered2srf(tnrb1, tnrb2, ed, pts1, pts2);
% 
% plot3(pnts1(1), pnts1(2), pnts1(3), 'ko');
% plot3(pnts2(1), pnts2(2), pnts2(3), 'k*');
% 
% plot3(tnrb1.points(ed,1), tnrb1.points(ed,2), tnrb1.points(ed,3), 'g.', 'MarkerSize', 15);
% plot3(tnrb1.points(ed,1), tnrb1.points(ed,2), tnrb1.points(ed,3), 'g', 'LineWidth', 1.5);
% 
% figure; hold on;
% triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2)); 
% plot(tnrb1.nodes(p1,1), tnrb1.nodes(p1,2), 'r*');
% plot(tnrb1.nodes(ed,1), tnrb1.nodes(ed,2), 'g.', 'MarkerSize', 15);
% plot(pts1(1), pts1(2), 'ro'); 
% title('Parametric mesh of surface 1'); 
% axis equal;
% 
% figure; hold on;
% triplot(tnrb2.delaunay, tnrb2.nodes(:,1), tnrb2.nodes(:,2)); 
% plot(tnrb2.nodes(p2,1), tnrb2.nodes(p2,2), 'r*');
% plot(pts2(1), pts2(2), 'ro'); 
% title('Parametric mesh of surface 2'); 
% axis equal;



