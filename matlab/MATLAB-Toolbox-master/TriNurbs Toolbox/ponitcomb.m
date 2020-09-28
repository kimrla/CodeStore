function [npts,npnts]=ponitcomb(pts1,pts2,pnts1,pnts2,sed,stri)
% Combine the 2 series of intersections obtained by tnrbinterconnct or tnrbintersort into 1
% series intersections.The input 2 series intersections have to be sorted
% firstly as the same continuous curve(eg. Both are anticlockwise).
%物理域上的交点结合后两个曲面是相同的水密的，但参数域上的交点不一样。交点所在的边集合要包含交点所在的三角形集合，即交点开始和结束处的三角形不能超过开始和结束处的边
% Input:
%   pts1: Parameter coordinates of intersections on the grids.
%   pts2: Parameter coordinates of intersections which is NOT on the grids
%       but inside the triangulars.
%   pnts1: Coordinates of intersections on the grids.
%   pnts2: Coordinates of intersections which is NOT on the grids
%       but inside the triangulars.
%   sed: Edges of intersections of tnrb1.
%   stri: Triangulars of intersections of tnrb1.
% Output:
%   npts: New parameter intersections combined by 2 series original points.
%   npnts: New coordinates of intersections combined by 2 series original points.


m=length(pnts1);
n=length(pnts2);
beginpt=0;
cp=sort(stri(1,:));
%寻找第一个交点所在的边,与第一个交点所在的三角形对应
for i=1:m-1
    pc=unique([sed(i,:),sed(i+1,:)]);
    if (cp==pc)
        beginpt=i;
        break;
    end   
end
if (beginpt==0)
    beginpt=m;
end

npts=pts1(beginpt,:);
npnts=pnts1(beginpt,:);
index=beginpt;
for i=2:n
    if (sort(stri(i,:))==sort(stri(i-1,:)))
        npts=[npts;pts2(i-1,:)];
        npnts=[npnts;pnts2(i-1,:)];
        %最后一个点如果单独处于一个三角形，则不计算在内
    else
        index=index-1;
        if (index<1)
            index=index+m;
        end
        npts=[npts;pts1(index,:);pts2(i,:);];
        npnts=[npnts;pnts1(index,:);pnts2(i,:);];
    end
end

% nn=length(npts);
% for i=2:nn
%     if (npts(i,:)==npts(i-1,:))
%         npts(i,:)=[];
%         npnts(i,:)=[];
%     end
% end
end



% % The mesh seed length (h0)
% h0=1;
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
% [sed1_, stri2_, spts1_, spts2_, spnts1_, spnts2_]=tnrbintersects(tnrb2, tnrb1, p2t2, p2t1);
% 
% 
% 
% 
% 
% % Connect and extend the intersections
% inter1=tnrbinterconnct(tnrb1, tnrb2, p2t1, p2t2, {sed1, stri2, spts1, spts2, spnts1, spnts2});
% inter2=tnrbinterconnct(tnrb2, tnrb1, p2t2, p2t1, {sed1_, stri2_, spts1_, spts2_, spnts1_, spnts2_});
% sed1=inter1{1}; stri2=inter1{2};
% spts1=inter1{3}; spts2=inter1{4};
% spnts1=inter1{5}; spnts2=inter1{6};
% % nr=length(sed1);
% sed1_=inter2{1}; stri2_=inter2{2};
% spts1_=inter2{3}; spts2_=inter2{4};
% spnts1_=inter2{5}; spnts2_=inter2{6};
% 
% 
% [npts,npnts]=ponitcomb(spts1{1},spts2_{1},spnts1{1},spnts2_{1},sed1{1},stri2_{1});
% 
% 
% % % Plot the results
% % figure; hold on; 
% % tnrbplot(tnrb1); 
% % tnrbplot(tnrb2); 
% % axis equal; view(3); 
% % title('Geometric grid'); 
% % for r=1:nr
% %     plot3(spnts1{r}(:,1), spnts1{r}(:,2), spnts1{r}(:,3), 'ro', 'LineWidth', 1); 
% %     plot3(spnts2{r}(:,1), spnts2{r}(:,2), spnts2{r}(:,3), 'y', 'LineWidth', 1); 
% % end
% % 
% % figure; hold on;
% % triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2));  
% % for r=1:nr
% %     plot(spts1{r}(:,1), spts1{r}(:,2), 'k.', 'MarkerSize', 13); 
% %     plot(spts1{r}(:,1), spts1{r}(:,2), 'r', 'LineWidth', 1); 
% % end
% % title('Parametric mesh of surface 1');  
% % axis equal; 
% % 
% figure; hold on; 
% triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2)); 
% plot(npts(:,1), npts(:,2), 'k.', 'MarkerSize', 13); 
% plot(npts(:,1), npts(:,2), 'r', 'LineWidth', 1); 
% 
% axis equal;












