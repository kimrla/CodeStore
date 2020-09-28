function [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri2tri(tnrb1, tnrb2, p2t2, ei1, trg1, trg2)

% tnrbintertri2tri: Get a new intersection through a known intersection
% 
% Calling Sequences:
% 
%       [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri2tri(tnrb1, tnrb2, p2t2, ei1, trg1, trg2)
% 
% INPUTS:
% 
%       tnrb1, tnrb2 - Triangular representation of nurbs surface.
% 
%       p2t2 - The relations from points to triangles of tri-nurbs 
%                  surface 2. See also tnrbpts2tri.
% 
%       ei1 - An edge of tnrb1.
%
%       trg1 - A triangle of tnrb1 that includes ei1.
%
%       trg2 - A triangle of tnrb2 that intersects with ei1.
%
% OUTPUT: 
%
%      ed1  -  Edges of the triangles of tri-nurbs surface 1 that
%                   intersected with tri-nurbs surface 2.
%
%      tri2   -  Triangles tri-nurbs surface 2 that intersected with 
%                  tri-nurbs surface 1.
%
%       pts1, pts2 - Parametric intersection points of the two surfaces.
%
%       pnts1, pnts2 - Intersection points of the two surfaces.
%

pe1=trg1(((trg1~=ei1(1)) & (trg1~=ei1(2))));
if length(pe1)>1
    error('The triangle does not include the edge!');
end
et1=sort([pe1, ei1(1); pe1, ei1(2)], 2);
ed1=[]; tri2=[]; pts1=[]; pts2=[]; pnts1=[]; pnts2=[]; 
t=0;
for i=1:2
    ed=et1(i,:);
    for j=1:3
        n=length(p2t2{trg2(j)});
        for k=1:n            
            tri=tnrb2.delaunay(p2t2{trg2(j)}(k),:); 
            [tpnts2, tpnts1, x]=tnrbinterline2tri(tnrb1, tnrb2, ed, tri);
            tr=tnrb2.nodes(tri, :);
            line=tnrb1.nodes(ed, :); 
            if ~isempty(x)
                tpts1x=(1-x(3,:))*line(1,1)+x(3,:)*line(2,1);
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
% h0=1.8;
% 
% % Create a nurbs sphere
% center=[5,5,4];
% circ=nrbcirc(4, center, 0, pi);
% srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% 
% % Transform a nurbs surface into triangular representation
% tnrb1=nrb2tri(srf1, h0); 
% tnrb2=nrb2tri(srf2, h0); 
% 
% % Get the edges of a tri-nurbs surface that intersected with another tri-nurbs surface
% p2t1=tnrbpts2tri(tnrb1);
% p2t2=tnrbpts2tri(tnrb2);
% 
% % Get the intersection points of two tri-nurbs surfaces and sort them
% [sed1, stri2, spts1, spts2, spnts1, spnts2]=tnrbintersects(tnrb1, tnrb2, p2t1, p2t2);
% nr=length(sed1);
% 
% % Get a new intersection through a known interection
% rr=3; k=1;
% ei1=sed1{rr}(k,:); 
% trg0=RemDuplicate([sed1{rr}(k,:), sed1{rr}(k+1,:)]')';
% tri=tnrbedge2tri(p2t1, ei1); 
% tri=tnrb1.delaunay(tri,:);
% d=DistanceMatrix(trg0, tri);
% trg1=tri(d~=0,:);
% trg2=stri2{rr}(k,:);
% [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri2tri(tnrb1, tnrb2, p2t2, ei1, trg1, trg2);
% 
% % Plot the results
% figure; hold on; 
% tnrbplot(tnrb1); 
% tnrbplot(tnrb2); 
% axis equal; view(3); 
% title('Geometric grid'); 
% for r=1:nr
%     plot3(spnts1{r}(:,1), spnts1{r}(:,2), spnts1{r}(:,3), 'ro', 'LineWidth', 1); 
%     plot3(spnts2{r}(:,1), spnts2{r}(:,2), spnts2{r}(:,3), 'y', 'LineWidth', 1); 
% end
% trg2=[trg2,trg2(1)];
% plot3(tnrb1.points(ei1,1), tnrb1.points(ei1,2), tnrb1.points(ei1,3), 'g', 'LineWidth', 2);
% plot3(tnrb2.points(trg2,1), tnrb2.points(trg2,2), tnrb2.points(trg2,3), 'r', 'LineWidth', 2);
% plot3(pnts1(:,1), pnts1(:,2), pnts1(:,3), 'ro'); 
% plot3(pnts2(:,1), pnts2(:,2), pnts2(:,3), 'r*'); 
% 
% figure; hold on;
% triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2));  
% for r=1:nr
%     plot(spts1{r}(:,1), spts1{r}(:,2), 'k.', 'MarkerSize', 13); 
%     plot(spts1{r}(:,1), spts1{r}(:,2), 'r', 'LineWidth', 1); 
% end
% trg1=[trg1,trg1(1)];
% plot(tnrb1.nodes(trg1,1), tnrb1.nodes(trg1,2), 'r', 'LineWidth', 2);
% plot(tnrb1.nodes(ei1,1), tnrb1.nodes(ei1,2), 'g', 'LineWidth', 2); 
% plot(pts1(:,1), pts1(:,2), 'ro', 'MarkerSize', 13); 
% title('Parametric mesh of surface 1');  
% axis equal; 
% 
% figure; hold on; 
% triplot(tnrb2.delaunay, tnrb2.nodes(:,1), tnrb2.nodes(:,2)); 
% for r=1:nr
%     plot(spts2{r}(:,1), spts2{r}(:,2), 'k.', 'MarkerSize', 13); 
%     plot(spts2{r}(:,1), spts2{r}(:,2), 'r', 'LineWidth', 1); 
% end
% plot(tnrb2.nodes(trg2,1), tnrb2.nodes(trg2,2), 'g', 'LineWidth', 2);
% plot(pts2(:,1), pts2(:,2), 'ro', 'MarkerSize', 13); 
% title('Parametric mesh of surface 2'); 
% axis equal;






