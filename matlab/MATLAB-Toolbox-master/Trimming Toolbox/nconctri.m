function newtri=nconctri(tri,pts,bndid)
% Obtain the triangulation in the available domain after removing triangulars outside the
% boundary. This function is only available for parameter domain.The order
% of calling functions is ndelaunaytri.m-npointinparsrf.m-transindex.m-ncontri.m

% Input:
%   tri: Triangulation connectivity list in or out of the available domain.
%       The indices are corresponding to just pts. See triangulation.ConnectivityList
%   pts: Coordinates of points in the available parameter domain. See
%       triangulation.Points{1}.Although the points are all in the
%       available domain, the triangulars may be out of the domain.
%   bndid: Indices of points in the boundary corresponding to pts (not nodes),see boundptsid.The first element is the same
%       as the last one.If the surface contains an inner holl,there's
%       [NAN,NAN] in the boundary and its index is 0.And the outer loop is
%       in anticlockwise order and the innter loop clockwise order.
% Output:
%   newtri: New triangulation structure, including connectivity list newtri.delaunay and
%       the coordinates of points newtri.points.

newtri.points=pts;
newtri.delaunay=tri;
%,relations between each point to triangulars, see tnrbpts2tri
nt=length(tri);
nd=length(pts);
p2t=cell(1, nd);
for ii=1:nt
    for jj=1:3
        p2t{tri(ii,jj)}=[p2t{tri(ii,jj)}, ii];
    end
end

tem1=find(bndid==0);
if isempty(tem1)% The original surface is trimmed into 2 surfaces.
    bndpts=pts(bndid,:);
    bndpts_=bndpts;
    bndidid=bndid;
else
    % if there is a inner loop, add [nan,nan] to the coordinates of points
    % in the inner and outer loop for calling the function inpolygon.m
    tem2=pts(bndid(1),:);
    if (tem2(1)~=0 && tem2(2)~=0 && tem2(1)~=1 && tem2(2)~=1)
        bndidid=bndid(1:tem1-1);
        bndpts=pts(bndidid,:);     
        bndpts_=[pts(bndid(tem1+1:end),:);nan,nan;bndpts];
    else
        bndidid=bndid(tem1+1:end);
        bndpts=pts(bndidid,:);      
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
        % whether the edge connecting from i to j should be removed
        if (sum(de)==1)                 
            even=(pts(bndidid(i),:)+pts(bndidid(j),:))/2;
            in=inpolygon(even(1),even(2),bndpts_(:,1)',bndpts_(:,2)');
            if (~in)% The midpoint is out of the boundary
                triid=[triid,tnrbedge2tri(p2t,[bndidid(i),bndidid(j)])];   
            end
        end
    end
end

cfd(triid)=0;
cfd=logical(cfd);
newtri.delaunay=newtri.delaunay(cfd,:);

end

%% demo
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
% tritol=h0/4;
% [npts_,npnts_,ntri,intid]=ndelaunaytri(tnrb1,npts,npnts2,tritol);
% 
% indexpts1=[intid(:,1);intid(end,2)]';
% [pts,in,on,boundpts,boundnodesid,boundptsid]=npointinparsrf(npts_,indexpts1);
% 
% nnpnts=cell(1,2);
% nnpnts{1}=npnts_(in,:);
% % nnpnts{2}=[npnts_(indexpts1,:);npnts_(~in,:)]; % 没有挖孔，把交点放在最前面
% nnpnts{2}=npnts_(~in | on,:); % 有挖孔
% 
% k=2; 
% in_=~in;
% in_(indexpts1)=1; 
% 
% % k=1;
% % in_=in;
% 
% newtri=transindex(ntri,in_); % newtri可能在有效域外还有三角形
% newtri_=nconctri(newtri,pts{k},boundptsid{k});
% 
% figure;hold on;
% triplot(newtri_.delaunay,pts{k}(:,1),pts{k}(:,2));
% % plot(boundpts{k}(:,1),boundpts{k}(:,2),'r--','Marker','*','MarkerSize',4);
% axis equal;
% 
% figure;hold on;
% triplot(newtri,pts{k}(:,1),pts{k}(:,2));
% % plot(boundpts{k}(:,1),boundpts{k}(:,2),'r--','Marker','*','MarkerSize',4);
% axis equal;
% 
% figure;hold on;
% trisurf(newtri_.delaunay,nnpnts{k}(:,1),nnpnts{k}(:,2),nnpnts{k}(:,3));
% axis equal;
% 
% figure;hold on;
% trisurf(newtri,nnpnts{k}(:,1),nnpnts{k}(:,2),nnpnts{k}(:,3));
% axis equal;




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    