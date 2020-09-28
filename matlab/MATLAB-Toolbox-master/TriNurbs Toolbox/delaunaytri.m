function [npts,npnts,ntri,intid]=delaunaytri(nrb,pts1,pts2,pnts1,pnts2,tritol)
% Mesh the parameter and physical domain after getting the intersections.
% The intersection curve must sepearate the original 1 surface into 2
% different surfaces.
%delaunayTriangulation会自己在边界上添加结点
% Input:
%   nrb: Nurbs surface structure.
%   pts1:Parameter coordinates of knots in the parameter domain of the
%       original surface
%   pts2:Parameter coordinates of intersection points.The first and last
%       points may not be coincident although they are the same ones in
%       pyhsical domain.
%   pnts1:Physical coordinates of nodes in the physical domain of the
%       original surface
%   pnts2:Physical coordinates of intersection points.The first and last
%       points are coincident so that it is a closed loop (intersection).
%   disttol:Tolearnce of distance between the first point and the last
%       point of the intersections, which is used to determine whether the
%       intersection curve is a continuous inner loop.
%   tritol:Tolerance of the edge of the triangulars.If the length of the
%       edge is smaller than tritol, this triangular need to be removed.
% Output:
%   npts:New parameter coordinates of all the points both in the original
%       surface and intersections
%   npnts:New physical coordinates of all the points both in the original
%       surface and intersections
%   ntri:New triangulation of all the points,whose triangulars are regular
%   intid:Indices of points in the intersections corresponding to all the
%       ponits npts or npnts

%Triangulation structure(tnrb1, tnrb2) includes the parameter and pyhsical
%coordinates and the indice of truanglars after triangulation of the
%trimming surface.The triangulation is based on delaunayTriangulation, but
%with removing 3 kinds of abnormal triangulars.After dealing the odd
%triangulars, the delaunayTriangulation will be processed again.
    
n=length(pts2);
npts=[pts2;pts1];
npnts=[pnts2;pnts1];
% pnts=unique(pnts,'rows','stable');
%不必消除重合点，因为默认交点不在参数域节点上
tempind=[];
bndin=[1:n-1;2:n]';

% The first and last points are the same in the parameter domain, which
% means the intersections are circled as a loop
%最终输出参数域和物理域坐标时，如果是内部有孔的情况，则孔边界的第一个点和最后一个点已经去掉不再重复
if (norm(pts2(1,:)-pts2(end,:))<eps)
    npts(n,:)=[];
    npnts(n,:)=[];
    bndin=[1:n-1;2:n-1,1]';
end

    
%如果交点围成一个圈，则将最后一个点和第一个点连接起来：需要修改，物理域是一个圈，但参数域不一定
%如果是点重合，则delaunayTriangulation划分时会自动删掉多余重合点
%不加边界限制（即没有第3个参数时），delaunayTriangulation会自动删除重复的点;
%加边界限制，即输入第3个参数时，下标索引要围成一个圈，且注意对应好，防止因自动删除点对应错误情况
% cishu=1;
while (true)
    tri=delaunayTriangulation(npts,bndin);%该函数可能会自己在边界上添加结点，导致前后结点不一样
    %Warning: Intersecting edge constraints have been split, this may have added new points into the triangulation. 
    %delaunayTriangulation自动利用或插入形成约束边界的节点，导致剖分的参数域和物理域结点不匹配
%     cishu=cishu+1;
    %determine whether the funtion delaunayTriangulation inserts points
    %into the boundary, or deduces new points into the boundary
    %进行新插入节点的处理时，默认交线非自交
    npts=tri.Points;
    tricnt=tri.ConnectivityList;
    bndin_=tri.Constraints;
    
    if (length(bndin)==length(bndin_))
        bndin=bndin_;
    else        
        [npnts,bndin]=delaunaydeal(nrb,npts,bndin_,npnts,bndin(1,:));

    end
    pand=unique(bndin);   
    numtri=length(tricnt);
    numpts=length(npts);
    for i=1:numtri
        l1=norm(npnts(tricnt(i,1),:)-npnts(tricnt(i,2),:));
        l2=norm(npnts(tricnt(i,3),:)-npnts(tricnt(i,2),:));
        l3=norm(npnts(tricnt(i,1),:)-npnts(tricnt(i,3),:));
        L=sort([l1,l2,l3]);
        %情况1，一条边极小的情况
        if (L(1)<tritol)
            if (l1==L(1))
                if (~sum(pand==tricnt(i,1)))%如果有新的边界下标，判定条件全部需要改过来，边界下标可能超过n，以下同理。后面的判定条件也许要改。
                    tempind=[tempind,tricnt(i,1)];
                else
                    tempind=[tempind,tricnt(i,2)];
                end
            elseif (l2==L(1))
                if (~sum(pand==tricnt(i,2)))
                    tempind=[tempind,tricnt(i,2)];
                else
                    tempind=[tempind,tricnt(i,3)];
                end
            else
                if (~sum(pand==tricnt(i,1)))
                    tempind=[tempind,tricnt(i,1)];
                else
                    tempind=[tempind,tricnt(i,3)];
                end
            end
        else
            %第二种情况,两边之和接近最长边，且最长边处于交线上，删掉不在最长边上的点。
            if (abs(L(1)+L(2)-L(3))<tritol)
                tembndin=sort(bndin,2);                               
                if (l1==L(3))
                    tem9=sort([tricnt(i,1),tricnt(i,2)]);
                    tem8=sum(tembndin(:,1)==tem9(1) & tembndin(:,2)==tem9(2));
                    if (tem8)
                        tempind=[tempind,tricnt(i,3)];
                    end
                elseif (l2==L(3))
                    tem9=sort([tricnt(i,3),tricnt(i,2)]);
                    tem8=sum(tembndin(:,1)==tem9(1) & tembndin(:,2)==tem9(2));
                    if (tem8)
                        tempind=[tempind,tricnt(i,1)];
                    end
                else
                    tem9=sort([tricnt(i,3),tricnt(i,1)]);
                    tem8=sum(tembndin(:,1)==tem9(1) & tembndin(:,2)==tem9(2));
                    if (tem8)
                        tempind=[tempind,tricnt(i,2)];
                    end
                end                              
            end
        end
    end
    if (~isempty(tempind))
        logind=ones(1,numpts);
        logind(tempind)=0;
        logind=logical(logind);
        npts=npts(logind,:);
        npnts=npnts(logind,:);
        tempind=[];
        continue;
    else
        break;
    end
end

%  tem1=sort(tricnt);
tembndin=sort(bndin,2); 
for i=1:numtri
    l1=norm(npnts(tricnt(i,1),:)-npnts(tricnt(i,2),:));
    l2=norm(npnts(tricnt(i,3),:)-npnts(tricnt(i,2),:));
    l3=norm(npnts(tricnt(i,1),:)-npnts(tricnt(i,3),:));
    L=sort([l1,l2,l3]);
    indtri=[0,0,0];
    %第三种情况，两边之和接近最长边，且最长边不处于交线上
    if (L(1)>=tritol && abs(L(1)+L(2)-L(3))<tritol)              
        if (l1==L(3))
            indtri(1)=tricnt(i,1);
            indtri(2)=tricnt(i,2);
            indtri(3)=tricnt(i,3);
        elseif (l2==L(3))
            indtri(1)=tricnt(i,2);
            indtri(2)=tricnt(i,3);
            indtri(3)=tricnt(i,1);
        else
            indtri(1)=tricnt(i,3);
            indtri(2)=tricnt(i,1);
            indtri(3)=tricnt(i,2);
        end        
    %         tem2=sort(indtri(1),indtri(2));
        tem9=sort([indtri(1),indtri(2)]);
        tem8=sum(tembndin(:,1)==tem9(1) & tembndin(:,2)==tem9(2));
        if (~tem8)
            for j=1:numtri
                if (j~=i)
                    tem3=logical(tricnt(j,:)==indtri(1) | tricnt(j,:)==indtri(2));
                    if (sum(tem3)==2)
                        tem4=tricnt(j,:);
                        tem5=tem4(~tem3);
                        tricnt(i,:)=sort([indtri(1),indtri(3),tem5]);
                        tricnt(j,:)=sort([indtri(2),indtri(3),tem5]);
                        break;
                    end
                end
            end
        end
    end
end

ntri=tricnt;
intid=bndin;

end

% % The mesh seed length (h0)
% h0=2;
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
% % cc=find(spts1{1}(:,1)==spts1{1}(2,1));
% % if (length(cc)~=1)
% %     cc=cc(2);
% %     sed1{1}(cc:end,:)=[];
% %     stri2{1}(cc:end,:)=[];
% %     spts1{1}(cc:end,:)=[];
% %     spts2{1}(cc:end,:)=[];
% %     spnts1{1}(cc:end,:)=[];
% %     spnts2{1}(cc:end,:)=[];
% % end

% dist=h0;
% [npts,npnts1,npnts2,ndist]=equinter(spnts1{1},tnrb1,dist);
% 
% tritol=h0/10;
% [npts_,npnts_,ntri,intid]=delaunaytri(srf1,tnrb1.nodes,npts,tnrb1.points,npnts2,tritol);
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
% 
% figure;hold on;
% triplot(ntri,npts_(:,1),npts_(:,2));
% mm=length(intid);
% for i=1:mm
%     plot([npts_(intid(i,1),1),npts_(intid(i,2),1)],[npts_(intid(i,1),2),npts_(intid(i,2),2)],'ro');
%     plot([npts_(intid(i,1),1),npts_(intid(i,2),1)],[npts_(intid(i,1),2),npts_(intid(i,2),2)],'g');
% end
% axis equal;
% figure;hold on;
% trisurf(ntri,npnts_(:,1),npnts_(:,2),npnts_(:,3));
% mm=length(intid);
% for i=1:mm
%     plot3([npnts_(intid(i,1),1),npnts_(intid(i,2),1)],[npnts_(intid(i,1),2),npnts_(intid(i,2),2)],[npnts_(intid(i,1),3),npnts_(intid(i,2),3)],'ro');
%     plot3([npnts_(intid(i,1),1),npnts_(intid(i,2),1)],[npnts_(intid(i,1),2),npnts_(intid(i,2),2)],[npnts_(intid(i,1),3),npnts_(intid(i,2),3)],'g');
% end
% axis equal;
% 
% 
% 
% 
% 
% 





