function [ed1, tri2, pts1, pts2, pnts1, pnts2]=qnrbinter(qnrb1, qnrb2)

% qnrbinter: Get the intersection points of two quad-nurbs surfaces
% 
% Calling Sequences:
% 
%       [ed1, tri2, pts1, pts2, pnts1, pnts2]=qnrbinter(tnrb1, tnrb2)
% 
% INPUTS:
% 
%       qnrb1, qnrb2 - Triangular representation of nurbs surface.
%
% OUTPUT: 
%
%      ed1  -  Edges of the quadrangles of quad-nurbs surface 1 that
%                   intersecte with triangles of quad-nurbs surface 2.
%
%      tri2   -  Triangles of quad-nurbs surface 2 that intersecte with 
%                  edges of quanrangles of quad-nurbs surface 1.
%
%       pts1, pts2 - Parametric intersection points of the two surfaces.
%
%       pnts1, pnts2 - Intersection points of the two surfaces.
%

% The nearest points of the two surfaces
tol=max([qnrb1.seeds(1), qnrb2.seeds(1)]);
[p1, p2]=nearpnts(qnrb1.UniqPt, qnrb2.UniqPt, tol);

% Get approximated intersections by quadrangles and triangles
nd=length(p1); t=0;
for k=1:nd
    qeds=uniqpt2qedges(qnrb1, p1(k));
    Tr=uniqpt2tri(qnrb2, p2(k));
    m=size(qeds, 1);
    n=length(Tr);
    for i=1:m
        for j=1:n
            tri=qnrb2.tri(Tr(j),:); 
            [tpnts2, tpnts1, x]=interline2intri(qnrb2.points(tri,:), qnrb1.points(qeds(i,:),:)); 
            tr=qnrb2.nodes(tri, :); 
            line=qnrb1.nodes(qeds(i,:), :); 
            if ~isempty(x)
                tpts1x=(1-x(3,:))*line(1,1)+x(3,:)*line(2,1);
                tpts1y=(1-x(3,:))*line(1,2)+x(3,:)*line(2,2);
                tpts1=[tpts1x(:), tpts1y(:)];
                tpts2=tripoint(tr(1,:), tr(2,:), tr(3,:), x(1,:), x(2,:));
                tk=size(x, 2);
                ed1(t+1:t+tk,:)=qeds(i,:);
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
[ed1, id]=RemDuplicate(sort(ed1, 2));
tri2=tri2(id,:);
pts1=pts1(id,:);
pts2=pts2(id,:);
pnts1=pnts1(id,:);
pnts2=pnts2(id,:);


%% demo
% % The mesh seed length (h0)
% h0=1.8;
% 
% % Create a nurbs sphere
% L=[7,10,0; 5,1,10];
% center=[5,5,4];
% circ=nrbcirc(4, center, 0, pi);
% srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% 
% % Transform a nurbs surface into quadrangular representation
% qnrb1=nrb2quad(srf1, h0);
% qnrb2=nrb2quad(srf2, h0);
% 
% % Get approximated intersections by quadrangles and triangles
% [ed1, tri2, pts1, pts2, pnts1, pnts2]=qnrbinter(qnrb1, qnrb2);
% 
% % Plot the surface and related quadrangle and line
% figure; hold on;
% trisurf(qnrb1.quad, qnrb1.points(:,1), qnrb1.points(:,2), qnrb1.points(:,3));
% trisurf(qnrb2.tri, qnrb2.points(:,1), qnrb2.points(:,2), qnrb2.points(:,3));
% plot3(pnts1(:,1), pnts1(:,2), pnts1(:,3), 'r*');
% plot3(pnts2(:,1), pnts2(:,2), pnts2(:,3), 'ro');
% view(3); axis equal;
% 
% % Plot the parameric domain of surface 1
% figure; hold on;
% quadplot(qnrb1, 'r');
% plot(pts1(:,1), pts1(:,2), 'r*');
% axis equal;
% 
% % Plot the parameric domain of surface 2
% figure; hold on;
% triplot(qnrb2.tri, qnrb2.nodes(:,1), qnrb2.nodes(:,2), 'r');
% plot(pts2(:,1), pts2(:,2), 'r*');
% axis equal;






