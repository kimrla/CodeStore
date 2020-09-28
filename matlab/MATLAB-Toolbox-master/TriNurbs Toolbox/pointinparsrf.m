function [pts,in,on,boundpts,boundptsid]=pointinparsrf(nodes,pts1,tol)
% Determine the points of each surface in parameter domain after trimming, sorting and removing process.

% Input:
%   nodes: parameter coordinates of all the points the surface, including points in the intersection curve(pts1), see tsrf.nodes.
%   pts1: parameter coordinates of the intersection curve, which have to be
%       sorted firstly.
%   tol: tolerance to determine whether the series of intersection points
%       can be connected into a holl.
% Output:
%   pts: parameter coordinates of the points in each surface which is a cell
%       array. 
%   in: indices of points that is in and on the boundary of trimming surface 1
%       parameter domain corresponding to all the nodes.
%   on: indices of points that is on the boundary of trimming surface 1
%       parameter domain corresponding to all the nodes.
%   boundpts: parameter coordinates of the points in boundary of the surface
%       or surfaces
%   boundptsid: index of the points on boundaries corresponding to all the
%       points in nodes

%   If length(pts)==1, the original surface is NOT trimmed into 2 different
%   sufaces and pts{1} stores the nodes and intersections of original surface.
%   If length(pts)==2, the original surface is trimmed into 2 surface and
%   pts{1},pts{2} stores the parameter coordinates of nodes in each
%   surface. Corresponding to the length of pts, bndpts{1} stores points in
%   the 4 edges of original surface and the intersection cuvre, otherwise,
%   bndpts{1},bndpts{2} store boundary points of the 2 trimming surfaces
%   respectively.

if (nargin==2)
    tol=eps;
end
n1=pts1(1,:);
n2=pts1(end,:);
uu(1)=n1(1);vv(1)=n1(2);
uu(2)=n2(1);vv(2)=n2(2);
[u,uid]=sort([uu(1),uu(2)]);
[v,vid]=sort([vv(1),vv(2)]);
fixu=fix(u);
fixv=fix(v);
n=[n1,n2];
m=sum(logical((abs(n-1)<eps | n<eps)));


nodesu1=nodes((nodes(:,1)==0),:);
[~,temuv]=sort(nodesu1(:,2),'descend');
nodesu1=nodesu1(temuv,:);
nodesu2=nodes((nodes(:,1)==1),:);
[~,temuv]=sort(nodesu2(:,2));
nodesu2=nodesu2(temuv,:);
nodesv1=nodes((nodes(:,2)==0),:);
[~,temuv]=sort(nodesv1(:,1));
nodesv1=nodesv1(temuv,:);
nodesv2=nodes((nodes(:,2)==1),:);
[~,temuv]=sort(nodesv2(:,1),'descend');
nodesv2=nodesv2(temuv,:);
nodesuv=[nodesv1;nodesu2;nodesv2;nodesu1];
nodesuv_=unique(nodesuv,'rows','stable');


%首末交点不能在曲面网格交点上
% The original is trimmed into 2 surfaces without holl and there's 16
% situations of the domain.
if (m==2)
    if ((abs(fixu(1)-u(1))<eps && abs(fixu(2)-u(2))<eps) || (abs(fixv(1)-v(1))<eps && abs(fixv(2)-v(2))<eps))
        if((abs(fixu(1)-u(1))<eps && abs(fixu(2)-u(2))<eps))
            if (abs(u(1)-u(2))<eps)%第一种情况
                if (abs(u(1))<eps)
                    bndpts{1}=nodesu1((nodesu1(:,2)<v(2) & nodesu1(:,2)>v(1)),:);
                    bndpts{2}=[nodesu1((nodesu1(:,2)<v(1)),:);nodesv1;nodesu2;nodesv2;nodesu1((nodesu1(:,2)>v(2)),:)];
                else
                    bndpts{1}=nodesu2((nodesu2(:,2)<v(2) & nodesu2(:,2)>v(1)),:);
                    bndpts{2}=[nodesu2((nodesu2(:,2)>v(2)),:);nodesv2;nodesu1;nodesv1;nodesu2((nodesu2(:,2)<v(1)),:)];
                end

            else%第二种情况
                bndpts{1}=[nodesu1(nodesu1(:,2)<vv(uid(1)),:);nodesv1;nodesu2(nodesu2(:,2)<vv(uid(2)),:)];
                bndpts{2}=[nodesu2(nodesu2(:,2)>vv(uid(2)),:);nodesv2;nodesu1(nodesu1(:,2)>vv(uid(1)),:)];

            end
        else
            if (abs(v(1)-v(2))<eps)%第三种情况
                if (abs(v(1))<eps)
                    bndpts{1}=nodesv1((nodesv1(:,1)<u(2) & nodesv1(:,1)>u(1)),:);
                    bndpts{2}=[nodesv1((nodesv1(:,1)>u(2)),:);nodesu2;nodesv2;nodesu1;nodesv1((nodesv1(:,2)<v(1)),:)];
                else
                    bndpts{1}=nodesv2((nodesv2(:,1)<u(2) & nodesv2(:,1)>u(1)),:);
                    bndpts{2}=[nodesv2((nodesv2(:,1)<u(1)),:);nodesu1;nodesv1;nodesu2;nodesv2((nodesv2(:,1)>u(2)),:)];
                end

            else%第四种情况
                bndpts{1}=[nodesv2(nodesv2(:,1)<uu(vid(2)),:);nodesu1;nodesv1(nodesv1(:,1)<uu(vid(1)),:)];
                bndpts{2}=[nodesv1(nodesv1(:,1)>uu(vid(1)),:);nodesu2;nodesv2(nodesv2(:,1)>uu(vid(2)),:)];

            end
        end
    else
        if (abs(u(1))<eps)
            if (abs(v(1))<eps)%第五种情况
                bndpts{1}=[nodesu1(nodesu1(:,2)<vv(uid(1)),:);nodesv1(nodesv1(:,1)<u(2),:)];
                bndpts{2}=[nodesv1(nodesv1(:,1)>u(2),:);nodesu2;nodesv2;nodesu1(nodesu1(:,2)>vv(uid(1)),:)];

            else%第六种情况
                bndpts{1}=[nodesv2(nodesv2(:,1)<u(2),:);nodesu1(nodesu1(:,2)>vv(uid(1)),:)];
                bndpts{2}=[nodesu1(nodesu1(:,2)<vv(uid(1)),:);nodesv1;nodesu2;nodesv2(nodesv2(:,1)>u(2),:)];

            end
        else
            if (abs(v(1))<eps)%第七种情况
                bndpts{1}=[nodesv1(nodesv1(:,1)>u(1),:);nodesu2(nodesu2(:,2)<vv(uid(2)),:)];
                bndpts{2}=[nodesu2(nodesu2(:,2)>vv(uid(2)),:);nodesv2;nodesu1;nodesv1(nodesv1(:,1)<u(1),:)];

            else%第八种情况
                bndpts{1}=[nodesu2(nodesu2(:,2)>v(1),:);nodesv2(nodesv2(:,1)>u(1),:)];
                bndpts{2}=[nodesv2(nodesv2(:,1)<u(1),:);nodesu1;nodesv1;nodesu2(nodesu2(:,2)<v(1),:)];

            end
        end

    end
    bndpts{1}=unique(bndpts{1},'rows','stable');
    bndpts{2}=unique(bndpts{2},'rows','stable');
    

    for iii=1:2
        temp1=norm(bndpts{iii}(1,:)-pts1(1,:));temp2=norm(bndpts{iii}(1,:)-pts1(end,:));
        if (temp1>temp2)
            boundpts{iii}=[pts1;bndpts{iii};pts1(1,:)];
        else
            boundpts{iii}=[pts1(end:-1:1,:);bndpts{iii};pts1(end,:)];
        end
    end
    
    [in,on]=inpolygon(nodes(:,1),nodes(:,2),boundpts{1}(:,1),boundpts{1}(:,2));
    %MATLAB in=inpolygon求出来的位置既包含在边界内的，也包含在边界上的
    pts{1}=nodes(in,:);
    pts{2}=[nodes(~in,:);pts1];
    
else
  
    % The original surface is not trimmed into 2 surfaces, or there is a holl as the inner boundary
    if (m==1)
        % The original surface is not trimmed into 2 surface, there's just an intersection curve in the surface       
        boundpts{1}=[nodesuv_;pts1];
        pts{1}=nodes;
        in=true(length(nodes),1);
        on=false(length(nodes),1);
        numbnd=length(boundpts{1});
        for ik=1:numbnd
            on=on | logical(pts{1}(:,1)==boundpts{1}(ik,1) & pts{1}(:,2)==boundpts{1}(ik,2));
        end
                
    else
        if (norm(n1-n2)<tol)
            %There's a holl in the surface
            boundpts{1}=pts1;%inner loop,the fisrt point is the same as the last one           
            boundpts{2}=[nodesuv_;nodesuv_(1,:);nan,nan;flipud(boundpts{1})];
            %outer surface boundary, which consists of both inner and outer loop, notice there is [nan,nan] in boundpts{2} of the boundary of srf2
            [in,on]=inpolygon(nodes(:,1),nodes(:,2),boundpts{1}(:,1),boundpts{1}(:,2));
            %in和on没问题，都是相对内环边界而言的
            pts{1}=nodes(in | on,:);
            pts{2}=nodes(~in | on,:);           
        else            
            %There's just an intersection curve in the surface and the curve is in the inner of parameter domain
            boundpts{1}=[nodesuv_;pts1];
            pts{1}=nodes; 
            in=true(length(nodes),1);
            on=false(length(nodes),1);
            numbnd=length(boundpts{1});
            for ik=1:numbnd
                on=on | logical(pts{1}(:,1)==boundpts{1}(ik,1) & pts{1}(:,2)==boundpts{1}(ik,2));
            end            
        end
    end
end
    
if (nargout==5)
    nn=length(boundpts);
    for ii=1:nn
        mm=length(boundpts{ii});
        for jj=1:mm
            if (~isnan(boundpts{ii}(jj,1)))    %if there is [nan,nan] in boundpts, the index in boundptsid corresponding to [nan,nan] is 0    
                temp=find(pts{ii}(:,1)==boundpts{ii}(jj,1) & pts{ii}(:,2)==boundpts{ii}(jj,2));%boundpts中第一个点和最后一个点重合
                boundptsid{ii}(jj)=temp;
            end
        end
    end
end
end
             
      
    
% h0=1;
% 
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
% %输出都是第一个点和最后一个点重复
% 
% dist=h0;
% [npts,npnts1,npnts2,ndist]=equinter(spnts1{1},tnrb1,dist);
% 
% tritol=h0/10;
% [npts_,npnts_,ntri,intid]=delaunaytri(srf1,tnrb1.nodes,npts,tnrb1.points,npnts2,tritol);%对应关系正常
%  %ntri没用到，重新划分网格多此一举，跟mibrmconctri里面的重新划分网格重复
% tempt=[intid(:,1);intid(end,2)];
% pts1=npts_(tempt,:);
% 
% [srfpts,in,on,boundpts,boundptsid]=pointinparsrf(npts_,pts1);
% 
% nnpnts=cell(1,2);
% nnpnts{1}=npnts_(in,:);
% nnpnts{2}=[npnts_(~in,:);npnts_(tempt,:)];
% 
% cp=1;
% %本实例中中点插入恢复边界法没发挥作用（即没有插入新的边界点），如果有新的边界点，必须在物理域上对应位置求出新的物理边界点
% tr=mibrmconctri(srfpts{cp},boundpts{cp},boundptsid{cp},nnpnts{cp});
% newtri=conctri(tr.delaunay,tr.nodes,boundptsid{cp});
% 
% figure;
% triplot(newtri.delaunay,newtri.points(:,1),newtri.points(:,2));
% axis equal;view(2);
%   
% figure;
% trisurf(newtri.delaunay,nnpnts{cp}(:,1),nnpnts{cp}(:,2),nnpnts{cp}(:,3));
% axis equal;view(2);
