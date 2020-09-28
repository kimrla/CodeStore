function newtri=conctri(tri,pts,bndid)
% Obtain the triangulation after removing triangulars outside the
% boundary.This function is used uaually after executing the function
% mibrmconctri.m. This function is only available for parameter domain.

% Input:
%   tri: Triangulation connectivity list. See triangulation.ConnectivityList
%   pts: Coordinates of points in parameter domain. See triangulation.Points
%   bndid: Indices of points in the boundary.The first element is the same
%       as the last one.If the surface contains an inner holl,there's
%       [NAN,NAN] in the boundary and its index is 0.
% Output:
%   newtri: New triangulation structure, including connectivity list and
%       the coordinates of points


newtri.delaunay=tri;
newtri.points=pts;

%,relations between points to triangulars, see tnrbpts2tri
nt=length(tri);%numbers of triangulars
nd=length(pts);%numbers of points
p2t=cell(1, nd);
for ii=1:nt
    for jj=1:3
        p2t{tri(ii,jj)}=[p2t{tri(ii,jj)}, ii];
    end
end

tem1=find(bndid==0);
if isempty(tem1)    
    bndpts=pts(bndid,:);%coordinates of points in boundary
    bndpts_=bndpts;
    bndidid=bndid;
else
    % if there is a inner loop, add [nan,nan] to the coordinates of points
    % in loop boundary
    tem2=pts(bndid(1),:);
    if (tem2(1)~=0 && tem2(2)~=0 && tem2(1)~=1 && tem2(2)~=1)
        bndidid=bndid(1:tem1-1);
        bndpts=pts(bndidid,:);     
        %inpolygon中，xv和yv内外环的转向要相反,默认外圈逆时针，内圈顺时针
        % if bndpts is anticlockwise
        panding1=bndpts(1,:)-bndpts(2,:);panding1(3)=0;
        panding2=bndpts(3,:)-bndpts(2,:);panding2(3)=0;
        panding3=cross(panding1,panding2);%顺正逆负
        if (sum(panding3)<0)
            bndpts=bndpts(end:-1:1,:);
        end
        bndpts_=[pts(bndid(tem1+1:end),:);nan,nan;bndpts];
    else
        bndidid=bndid(tem1+1:end);
        bndpts=pts(bndidid,:);      
        %inpolygon中，xv和yv内外环的转向要相反,默认外圈逆时针，内圈顺时针
        % if bndpts is anticlockwise
        panding1=bndpts(1,:)-bndpts(2,:);panding1(3)=0;
        panding2=bndpts(3,:)-bndpts(2,:);panding2(3)=0;
        panding3=cross(panding1,panding2);%顺正逆负
        if (sum(panding3)<0)
            bndpts=bndpts(end:-1:1,:);
        end
        bndpts_=[pts(bndid(1:tem1-1),:);nan,nan;bndpts];
    end
end

n=length(bndpts);
triid=[];
cfd=ones(1,length(tri));
for i=1:n-3
    for j=i+2:n-1
        %relations between points to edges, see tnrbpt2edges
        pt=bndidid(i);
        temtri=p2t{pt};
        nnn=length(temtri);
        edges=zeros(2*nnn, 1);
        for ii=1:nnn
            p=tri(temtri(ii),:)~=pt;
            edges(2*ii-1:2*ii)=tri(temtri(ii),p);
        end  
        edges=RemDuplicate(edges);
                       
        de=logical(edges==bndidid(j));
        %whether the edge connecting from i to j exists
        if (sum(de)==1)                 
            even=(pts(bndidid(i),:)+pts(bndidid(j),:))/2;
            in=inpolygon(even(1),even(2),bndpts_(:,1)',bndpts_(:,2)');
            if (~in)%中点落在边界外
                triid=[triid,tnrbedge2tri(p2t,[bndidid(i),bndidid(j)])];   
            end
        end
    end
end

cfd(triid)=0;
cfd=logical(cfd);
newtri.delaunay=newtri.delaunay(cfd,:);
                
%% demo1
% h0=0.5;
% 
% % Create a plane square and a plane cuve
% lin1=nrbline([0,1], [9,1]);
% lin2=nrbline([0,6], [9,6]);
% srf=nrbruled(lin1, lin2);
% crv=nrbtestcrv;
% 
% % Transform a nurbs surface into triangular representation
% tsrf=nrb2tri(srf, h0);
% tcrv=nrb2tri(crv, h0);
% 
% % The nearest points from the surface to the curve
% tol=max([tsrf.seeds(1), tcrv.seeds(1)]);
% [p1, p2, d]=nearpnts(tsrf.points, tcrv.points, tol);
% 
% % Get the relations from points to triangles of tri-nurbs
% p2t1=tnrbpts2tri(tsrf);
% 
% % Get the intersection points of a tri-nurbs surface with a curve
% [ed1, ed2, pts1, pts2, pnts1, pnts2]=tnrbinterline(tsrf, tcrv, p2t1, p1, p2);       
% 
% dist=h0;
% [npts,npnts1,npnts2,ndist]=equinter(pnts1,tsrf,dist);%处理任意曲面srf2
% 
% tritol=h0/10;
% [npts_,npnts_,~,intid]=delaunaytri(srf,tsrf.nodes,npts,tsrf.points,npnts2,tritol);%对应关系正常
% 
% %交点对应的下标，第1个和最后一个交点相同
% tempt=[intid(:,1);intid(end,2)];
% pts1=npts_(tempt,:);
% 
% [srfpts,in,on,boundpts,boundptsid]=pointinparsrf(npts_,pts1);
% 
% nnpnts=cell(1,2);
% nnpnts{1}=npnts_(in | on,:);
% nnpnts{2}=[npnts_(~in,:);npnts_(tempt,:)];
% 
% cp=2;
% tr=mibrmconctri(srfpts{cp},boundpts{cp},boundptsid{cp},nnpnts{cp});
% 
% newtri=conctri(tr.delaunay,tr.nodes,boundptsid{cp});
% 
% figure;
% triplot(newtri.delaunay,newtri.points(:,1),newtri.points(:,2));
% axis equal;view(2);
%   
% figure;
% trisurf(newtri.delaunay,tr.points(:,1),tr.points(:,2),tr.points(:,3));
% axis equal;view(2);

%% demo2
% clear;
% h0=1;
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
% % Connect and extend the intersections
% inter1=tnrbinterconnct(tnrb1, tnrb2, p2t1, p2t2, {sed1, stri2, spts1, spts2, spnts1, spnts2});
% sed1=inter1{1}; stri2=inter1{2};
% spts1=inter1{3}; spts2=inter1{4};
% spnts1=inter1{5}; spnts2=inter1{6};
% 
% dist=h0;
% [npts,npnts1,npnts2,ndist]=equinter(spnts2{1},tnrb2,dist);%处理任意曲面srf2
% 
% tritol=h0/10;
% [npts_,npnts_,~,intid]=delaunaytri(srf2,tnrb2.nodes,npts,tnrb2.points,npnts2,tritol);%对应关系正常
% 
% %交点对应的下标，第1个和最后一个交点相同
% tempt=[intid(:,1);intid(end,2)];
% pts1=npts_(tempt,:);
% 
% [srfpts,in,on,boundpts,boundptsid]=pointinparsrf(npts_,pts1);
% 
% nnpnts=cell(1,2);
% nnpnts{1}=npnts_(in | on,:);
% tempt(end)=[];
% nnpnts{2}=npnts_(~in | on,:);   % %nnpnts{2}=[npnts_(~in,:);npnts_(tempt,:)];非挖孔情况下的点的对应关系

% cp=2;
% tr=mibrmconctri(srfpts{cp},boundpts{cp},boundptsid{cp},nnpnts{cp});
% 
% newtri=conctri(tr.delaunay,tr.nodes,boundptsid{cp});
% 
% figure;
% triplot(newtri.delaunay,newtri.points(:,1),newtri.points(:,2));
% axis equal;view(2);
%   
% figure;
% trisurf(newtri.delaunay,tr.points(:,1),tr.points(:,2),tr.points(:,3));
% axis equal;view(2);

%% demo
% clear;
% h0=1;
% center=[5,0,4];
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
% % Connect and extend the intersections
% inter1=tnrbinterconnct(tnrb1, tnrb2, p2t1, p2t2, {sed1, stri2, spts1, spts2, spnts1, spnts2});
% sed1=inter1{1}; stri2=inter1{2};
% spts1=inter1{3}; spts2=inter1{4};
% spnts1=inter1{5}; spnts2=inter1{6};
% 
% dist=h0;
% [npts,npnts1,npnts2,ndist]=equinter(spnts1{1},tnrb1,dist);%处理任意曲面srf2
% 
% tritol=h0/10;
% 
% [npts_,npnts_,~,intid]=delaunaytri(srf1,tnrb1.nodes,npts,tnrb1.points,npnts2,tritol);%对应关系正常
% 
% %交点对应的下标，第1个和最后一个交点相同
% tempt=[intid(:,1);intid(end,2)];
% pts1=npts_(tempt,:);
% 
% [srfpts,in,on,boundpts,boundptsid]=pointinparsrf(npts_,pts1);
% 
% nnpnts=cell(1,2);
% nnpnts{1}=npnts_(in | on,:);
% nnpnts{2}=[npnts_(~in,:);npnts_(tempt,:)];
% 
% cp=2;
% tr=mibrmconctri(srfpts{cp},boundpts{cp},boundptsid{cp},nnpnts{cp});
% 
% newtri=conctri(tr.delaunay,tr.nodes,boundptsid{cp});
% 
% figure;
% triplot(newtri.delaunay,newtri.points(:,1),newtri.points(:,2));
% axis equal;view(2);
%   
% figure;
% trisurf(newtri.delaunay,tr.points(:,1),tr.points(:,2),tr.points(:,3));
% axis equal;view(2);


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    