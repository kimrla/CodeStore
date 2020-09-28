function [npts,npnts,ntri,intid]=ndelaunaytri(tnrb,pts2,pnts2,tritol)
% Mesh the parameter and physical domain after getting the intersections.
% The intersection curve must sepearate the original 1 surface into 2
% different surfaces.

% Input:
%   tnrb: Tri-nurbs surface structure.
%   pts2:Parameter coordinates of intersection points.The first and last
%       points may not be coincident although they are the same ones in
%       pyhsical domain.
%   pnts2:Physical coordinates of intersection points.The first and last
%       points are coincident so that it is a closed loop (intersection).
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

nrb=tnrb.nurbs;
pts1=tnrb.nodes;%pts1:Parameter coordinates of knots in the parameter domain of the tri-surface corresponding to the NURBS surface nrb.
pnts1=tnrb.points;%pnts1:Physical coordinates of nodes in the physical domain of the tri-surface corresponding to the NURBS surface nrb.

n=length(pts2);
npts=[pts2;pts1];
npnts=[pnts2;pnts1];

bndin=[1:n-1;2:n]';
% Considering the function delaunayTriangulation will delete the same points,if the first and last points are the same in the parameter domain, which
% means the intersections are circled as a loop,then remove the parametric coordinates of the last point in npts/npnts and modify bndin.
if (norm(pts2(1,:)-pts2(end,:))<eps)
    npts(n,:)=[];
    npnts(n,:)=[];
    bndin=[1:n-1;2:n-1,1]';% construct the closed boundary compulsorily
end

tempind=[];
while (true)
    % This function may insert or delete nodes automatically if the boundary condition (bndin) is not regular.
    tri=delaunayTriangulation(npts,bndin);
    %Determine whether the funtion delaunayTriangulation inserts points
    %into the boundary, or deduces new points into the boundary
    npts=tri.Points;
    tricnt=tri.ConnectivityList;
    bndin_=tri.Constraints;
    
    % This procedure is to deal with the condition that the function
    % delaunayTriangulation inserted new points or induced nodes which
    % pass the boundary.In that case, npts may consist of new rows 
    % (because of inserting new points) and the boundary indices bndin_ is not we want.
    % After the function delaunaydeal, the new boundary indices bndin whose
    % elements correspond to the real indices of boundary points in npts
    % and the new physical coordinates npnts corresponding to npts and bndin is obtained.
    if (length(bndin)==length(bndin_))
        bndin=bndin_;
    else        
        [npnts,bndin]=delaunaydeal(nrb,npts,bndin_,npnts,bndin);
    end
    
    pand=unique(bndin);   
    numtri=length(tricnt);
    for i=1:numtri
        l(1)=norm(npnts(tricnt(i,1),:)-npnts(tricnt(i,2),:));
        l(2)=norm(npnts(tricnt(i,3),:)-npnts(tricnt(i,2),:));
        l(3)=norm(npnts(tricnt(i,1),:)-npnts(tricnt(i,3),:));
        L=sort(l);
        
        % Situation 1: At least 1 edge is too small
        if (L(2)<tritol)
            tempind=[tempind,tricnt(i,1),tricnt(i,2)];% All 3 edges of this triangle are too small, then remove the first and second point.
        elseif (L(1)<tritol)
            for j=1:3
                if (l(j)==L(1))
                    if (~sum(pand==tricnt(i,j)))
                        tempind=[tempind,tricnt(i,j)];
                        break;
                    else
                        j=j+1;
                        if (j>3)
                            j=1;
                        end
                        if (~sum(pand==tricnt(i,j)))
                            tempind=[tempind,tricnt(i,j)];
                            break;
                        end
                    end
                end
            end                                       
        else
            % Situation 2: In one triangle,the sum of 2 edges is similar to the 3rd one and
            % the 3rd edge belongs to the intersection curve,then remove
            % the point which does not belongs to the 3rd edge.
            if (abs(L(1)+L(2)-L(3))<tritol)
                tembndin=sort(bndin,2);     
                for j=1:3
                    if (l(j)==L(3))
                        jj=j+1;
                        if (jj>3)
                            jj=1;
                        end
                        tem1=sort([tricnt(i,j),tricnt(i,jj)]);
                        tem2=sum(tembndin(:,1)==tem1(1) & tembndin(:,2)==tem1(2));
                        if (tem2)
                            jj=jj+1;
                            if (jj>3)
                                jj=1;
                            end
                            tempind=[tempind,tricnt(i,jj)];
                            break;
                        end
                    end
                end
            end
        end
    end
    if (~isempty(tempind))
        logind=ones(1,length(npts));
        logind(tempind)=0;
        logind=logical(logind);
        npts=npts(logind,:);
        npnts=npnts(logind,:);
        
        % Once some points in npts/npnts are removed, the indices of points
        % after the removing points in bndin have to be decreased.
        numtem=length(tempind);
        for i1=1:numtem
            detertem=tempind(i1);
            zzbone=bndin>detertem;
            bndin=bndin-zzbone;
        end

        tempind=[];
        continue;
    else
        break;
    end
end

tembndin=sort(bndin,2); 
for i=1:numtri
    l(1)=norm(npnts(tricnt(i,1),:)-npnts(tricnt(i,2),:));
    l(2)=norm(npnts(tricnt(i,3),:)-npnts(tricnt(i,2),:));
    l(3)=norm(npnts(tricnt(i,1),:)-npnts(tricnt(i,3),:));
    L=sort(l);
    indtri=[0,0,0];
    % Situation 3: In one triangle,the sum of 2 edges is similar to the 3rd one and
    % the 3rd edge doesn't belong to the intersection curve,then change the
    % diagonal of the quadrilateral corresponding to the triangle.
    if (abs(L(1)+L(2)-L(3))<tritol)  
        for k=1:3
            if (l(k)==L(3))
                k1=k+1;k2=k-1;
                if (k1>3)
                    k1=1;
                end
                if (k2<1)
                    k2=3;
                end
                indtri(1)=tricnt(i,k);
                indtri(2)=tricnt(i,k1);
                indtri(3)=tricnt(i,k2);
                break;
            end
        end
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
%% demo1
% % The mesh seed length (h0)
% h0=1.6;
% 
% % Create a nurbs sphere
% center=[5,5,4];
% circ=nrbcirc(4, center, 0, pi);
% srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
% srf2=nrbtestsrf;% change srf1 and srf2
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

% dist=h0;
% [npts,npnts1,npnts2,ndist]=nequinter(spnts1{1},tnrb1,dist);
% 
% tritol=h0/10;
% [npts_,npnts_,ntri,intid]=ndelaunaytri(tnrb1,npts,npnts2,tritol);
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
% plot3(spnts2{1}(:,1),spnts2{1}(:,2),spnts2{1}(:,3),'g*','MarkerSize',7);
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
% plot3(npnts2(:,1),npnts2(:,2),npnts2(:,3),'g*','MarkerSize',7);
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


%% demo2
% h0=1;
% 
% % Create a nurbs sphere
% center=[5,5,4];
% circ=nrbcirc(4, center, 0, pi);
% srf2=nrbrevolve(circ, center, [1,0,0], 2*pi);
% srf1=nrbtestsrf;
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
% cc=find(spts1{1}(:,1)==spts1{1}(2,1) & spts1{1}(:,2)==spts1{1}(2,2));
% if (length(cc)~=1)
%     cc=cc(2);
%     sed1{1}(cc:end,:)=[];
%     stri2{1}(cc:end,:)=[];
%     spts1{1}(cc:end,:)=[];
%     spts2{1}(cc:end,:)=[];
%     spnts1{1}(cc:end,:)=[];
%     spnts2{1}(cc:end,:)=[];
% end
% 
% dist=h0;
% [npts,npnts1,npnts2,ndist]=equinter(spnts1{1},tnrb1,dist);
% %%
% tritol=h0/4;
% [npts_,npnts_,ntri,intid]=ndelaunaytri(tnrb1,npts,npnts2,tritol);
% 
% %%
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
% plot3(spnts2{1}(:,1),spnts2{1}(:,2),spnts2{1}(:,3),'g*','MarkerSize',7);
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
% plot3(npnts2(:,1),npnts2(:,2),npnts2(:,3),'g*','MarkerSize',7);
% %%
% figure;hold on;
% triplot(ntri,npts_(:,1),npts_(:,2));
% mm=length(intid);
% for i=1:mm
%     plot([npts_(intid(i,1),1),npts_(intid(i,2),1)],[npts_(intid(i,1),2),npts_(intid(i,2),2)],'ro');
%     plot([npts_(intid(i,1),1),npts_(intid(i,2),1)],[npts_(intid(i,1),2),npts_(intid(i,2),2)],'g');
% end
% axis equal;
% 
% figure;hold on;
% trisurf(ntri,npnts_(:,1),npnts_(:,2),npnts_(:,3));
% mm=length(intid);
% for i=1:mm
%     plot3([npnts_(intid(i,1),1),npnts_(intid(i,2),1)],[npnts_(intid(i,1),2),npnts_(intid(i,2),2)],[npnts_(intid(i,1),3),npnts_(intid(i,2),3)],'ro');
%     plot3([npnts_(intid(i,1),1),npnts_(intid(i,2),1)],[npnts_(intid(i,1),2),npnts_(intid(i,2),2)],[npnts_(intid(i,1),3),npnts_(intid(i,2),3)],'g');
% end
% axis equal;

