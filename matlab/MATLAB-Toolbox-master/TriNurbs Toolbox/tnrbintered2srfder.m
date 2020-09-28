function [d, pts1, pnts1, pts2, pnts2]=tnrbintered2srfder(tnrb1, tnrb2, ed, pts1, pts2, it)

% tnrbintered2srfder: Get the intersection of an edge of a tri-nurbs surface to another tri-nurbs surfaces
% 
% Calling Sequences:
% 
%       [d, pts1, pnts1, pts2, pnts2]=tnrbintered2srfder(tnrb1, tnrb2, ed, pts1, pts2)
% 
%       [d, pts1, pnts1, pts2, pnts2]=tnrbintered2srfer(tnrb1, tnrb2, ed, pts1, pts2, it)
% 
% INPUTS:
% 
%       tnrb1, tnrb2 - Two tri-nurbs surfaces.
% 
%       ed -  Indexes of an edge on tnrb1.
% 
%       pts1, pts2 - Approximated parametric intersections of the two surfaces.
%
%       it -  The number of iterations.
%
% OUTPUT: 
%
%       d  -  Distances between the two sets intersection points.
%
%       pts1, pts2 - Parametric intersections of the two surfaces.
%
%       pnts1, pnts2 - Intersection points of the two surfaces.
%

if nargin==5
    it=3; 
end

edgseed = @(tnrb, r, ei, k) (1-r)*tnrb.nodes(ei(1), k)+r*tnrb.nodes(ei(2), k);
u1=tnrb1.nodes(ed(1), 1);
u2=tnrb1.nodes(ed(2), 1);
v1=tnrb1.nodes(ed(1), 2);
v2=tnrb1.nodes(ed(2), 2);
du=u2-u1; dv=v2-v1; 
if du~=0
    x=[pts2(:); (pts1(1)-u1)/du];
else
    x=[pts2(:); (pts1(2)-v1)/dv];
end
pp=x>1; x(pp)=1; 
pp=x<0; x(pp)=0; 
dnrb1=nrbderiv(tnrb1.nurbs);
dnrb2=nrbderiv(tnrb2.nurbs);
for i=1:it
    pts1(1)=edgseed(tnrb1, x(3), ed, 1);
    pts1(2)=edgseed(tnrb1, x(3), ed, 2);
    [pnts2, jac2]=nrbdeval(tnrb2.nurbs, dnrb2, x(1:2)); 
    [pnts1, jac1]=nrbdeval(tnrb1.nurbs, dnrb1, pts1(:));  
    jac1=du*jac1{1}+dv*jac1{2};
    dr=pnts2-pnts1; 
    F(1,1)=dot(dr, jac2{1}); 
    F(2,1)=dot(dr, jac2{2}); 
    F(3,1)=-dot(dr, jac1); 
    dF(1,1)=dot(jac2{1}, jac2{1}); 
    dF(1,2)=dot(jac2{1}, jac2{2}); 
    dF(1,3)=-dot(jac2{1}, jac1); 
    dF(2,1)=dF(1,2); dF(3,1)=dF(1,3); 
    dF(2,2)=dot(jac2{2}, jac2{2}); 
    dF(2,3)=-dot(jac2{2}, jac1); 
    dF(3,2)=dF(2,3); 
    dF(3,3)=dot(jac1, jac1); 
    x=x-dF\F; 
    pp=x>1; x(pp)=1; 
    pp=x<0; x(pp)=0; 
end
pts1(1)=edgseed(tnrb1, x(3), ed, 1);
pts1(2)=edgseed(tnrb1, x(3), ed, 2);
pts2=x(1:2);
pnts1=nrbeval(tnrb1.nurbs, pts1(:));  
pnts2=nrbeval(tnrb2.nurbs, pts2(:)); 
d=norm(pnts2-pnts1);


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
% k=58; i=2; 
% pt1=p1(k); pt2=p2(k); 
% p2t=tnrbpts2tri(tnrb1);
% edges=tnrbpt2edges(tnrb1, p2t, pt1);
% ed=[pt1, edges(i)];
% 
% % Get intersections of a edge of a tri-nurbs surface with another tri-nurbs surface
% pts1=(tnrb1.nodes(ed(1),:)+tnrb1.nodes(ed(2),:))/2;
% pts2=tnrb2.nodes(pt2, :);
% [dt, pts1, pnts1, pts2, pnts2]=tnrbintered2srfder(tnrb1, tnrb2, ed, pts1, pts2);
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





