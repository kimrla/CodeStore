function [npts,npnts1,npnts2,ndist]=nequinter(pnts,tnrb,dist)
%Get the paramater and physical coordinates of intersections that are
% equidistant.The interpolation uses NURL basis and then is transformed into
% NURBS structure.

% Input:
%   pnts:Physical coordinates of interpolation.
%   tnrb:Triangulation representation of a surface data structure.
%   dist:Distance between the equidistant intersections
% Output:
%   npts:New parameter coordinates of intersections of the surface that is
%       equidistantly distributed in the physical domain.
%   npnts1:New physical coordinates of intersections of the surface that is
%       equidistantly distributed, which are obtained by interpolation in
%       the NURL intersection curve.
%   npnts2:New physical coordinates of intersections of the surface that is
%       equidistantly distributed, which are obtained by reverse of the
%       NURBS surface, see nrbsrfreverse.
%   ndist:Distance between npnts1 and npnts2, which is used to check the
%       precision

% The watertight will be destroyed after calling on this function.

m=length(pnts);
l=[];
for i=1:m-1
    l=[l,norm(pnts(i,:)-pnts(i+1,:))];
end
L=sum(l);

pts1=0;
lsum=0;
for i=1:m-1
    lsum=lsum+l(i);
    pts1=[pts1,lsum/L];
end


nrl1=nrlmake(pnts',pts1);
n=fix(L/dist);
pts2=linspace(0,1,n+1);

npnts1=nrleval(nrl1,pts2);
npnts1=npnts1(:,:)';

srf=tnrb.nurbs;
nodes=tnrb.nodes;
points=tnrb.points;
[~,p2]=nearpnts(npnts1,points);

np=length(npnts1);
npnts2=zeros(np,3);
npts=zeros(np,2);
ndist=zeros(np,1);
dsrf=nrbderiv(srf);

for k=1:np
    [npts(k,:), npnts2(k,:), ndist(k)]=nrbsrfreverse(srf, dsrf, nodes(p2(k),:), npnts1(k,:)');
    % calculate the physical equidistant points in the NURBS surface more
    % explicitly (comparing to npnts1) and reverse the parameter coordinates of the points
end

%The fisrt point is the same as the last one
% if the first point is the same as the last one (which means the curve is
% a loop),then the parametric coordinates of reversing procedure of these 2
% points are the same,which needs to be composed.

if (norm(pnts(1,:)-pnts(end,:))<1e-5)
    det=logical(npts(1,:)<eps);
    if (det(1)==1)
        if (abs(npts(2,1))<abs(npts(2,1)-1))
            npts(end,1)=1;
        else
            npts(1,1)=1;
        end
    elseif (det(2)==1)
        if (abs(npts(2,2))<abs(npts(2,2)-1))
            npts(end,2)=1;
        else
            npts(1,2)=1;
        end
    end
end

end


% % The mesh seed length (h0)
% h0=1.5;
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
% 
% 
% % Connect and extend the intersections
% inter1=tnrbinterconnct(tnrb1, tnrb2, p2t1, p2t2, {sed1, stri2, spts1, spts2, spnts1, spnts2});
% sed1=inter1{1}; stri2=inter1{2};
% spts1=inter1{3}; spts2=inter1{4};
% spnts1=inter1{5}; spnts2=inter1{6};
% 
% figure; hold on; 
% triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2)); 
% axis equal;
% plot(spts1{1}(:,1), spts1{1}(:,2), 'k.', 'MarkerSize', 13); 
% plot(spts1{1}(:,1), spts1{1}(:,2), 'r', 'LineWidth', 1); 
% 
% figure;hold on;
% trisurf(tnrb1.delaunay,tnrb1.points(:,1),tnrb1.points(:,2),tnrb1.points(:,3));
% axis equal;
% plot3(spnts1{1}(:,1),spnts1{1}(:,2),spnts1{1}(:,3),'ro','MarkerSize',7);
% plot3(spnts1{1}(:,1),spnts1{1}(:,2),spnts1{1}(:,3),'r','LineWidth',1);
% plot3(spnts2{1}(:,1),spnts2{1}(:,2),spnts2{1}(:,3),'g*','MarkerSize',7);
% plot3(spnts1{1}(:,1),spnts1{1}(:,2),spnts1{1}(:,3),'g','LineWidth',1);
% 
% dist=h0/3;
% [npts,npnts1,npnts2,ndist]=equinter(spnts1{1},tnrb1,dist);
% 
% figure;hold on;
% triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2)); 
% axis equal;
% plot(npts(:,1), npts(:,2), 'k.', 'MarkerSize', 13); 
% plot(npts(:,1), npts(:,2), 'r', 'LineWidth', 1); 
% 
% figure;hold on;
% trisurf(tnrb1.delaunay,tnrb1.points(:,1),tnrb1.points(:,2),tnrb1.points(:,3));
% axis equal;
% plot3(npnts1(:,1),npnts1(:,2),npnts1(:,3),'ro','MarkerSize',7);
% plot3(npnts1(:,1),npnts1(:,2),npnts1(:,3),'r','LineWidth',1);
% plot3(npnts2(:,1),npnts2(:,2),npnts2(:,3),'g*','MarkerSize',7);
% plot3(npnts2(:,1),npnts2(:,2),npnts2(:,3),'g','LineWidth',1);



