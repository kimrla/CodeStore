function [pts,in,on,boundpts,boundnodesid,boundptsid]=npointinparsrf(nodes,indexpts1,tol)
% Determine the points of each surface in parameter domain after trimming,
% sorting and removing process. Notice that this function is only available
% to regular parameter domian.For degenerated domian,eg.par-domain of
% revolving sphere, this function should be modified.

% Input:
%   nodes: parameter coordinates of all the nodes the surface, including nodes in the intersection curve(pts1).
%   indexpts1: indices of the intersection curve, which needs be
%       sorted firstly.If the intersections are connected into a closed
%       loop, the fisrt element of pts1 is the same as the last one.The
%       intersection points pts1=nodes(indexpts1,:).
%   tol: tolerance to determine whether the series of intersection points
%       can be connected into a loop.
% Output:
%   pts: parameter coordinates of the points in each surface which is a cell
%       array. 
%   in: indices of points that is in/on the boundary of trimming surface 1
%       corresponding to all the nodes.Therefore, indices of points in
%       trimming surface 2 is ~in&pts1.
%   on: indices of points that is on the boundary of trimming surface 1
%       corresponding to all the nodes but they are NOT stored in a specific order.Therefore, indices of points only in
%       the domain of surface 1 is in&~on
%   boundpts: parameter coordinates of the points in boundaries of the
%       surfaces,which are sorted in a specific order,a cell array.If
%       there's a holl in the surface, there's [NAN,NAN] point in the
%       boundary to separate the inner and outer boundary.
%   boundnodesid: indices of points corresponding to nodes,which means
%       boundpts=nodes(boundnodesid,:),a cell array.If there's an inner holl,
%       the index corresponding to [NAN,NAN] is 0.
%   boundptsid: indices of points corresponding to pts,which means
%       boundpts=pts(boundptsid,:),a cell array.If there's an inner holl,
%       the index corresponding to [NAN,NAN] is 0.

%   If length(pts)==1, the original surface is NOT trimmed into 2 different
%   sufaces and pts{1} stores the nodes and intersections of original surface.
%   If length(pts)==2, the original surface is trimmed into 2 surface and
%   pts{1},pts{2} stores the parameter coordinates of nodes in each
%   surface. And acquiesce that all the intersection points are NOT
%   in coincidence with the nodes of surface's parameter domain.

if (nargin==2)
    tol=eps;
end

% Set the ending nodes of intersections and 4 boundary edges of the
% original surface's parameter domain.
pts1=nodes(indexpts1,:);
firstnode=pts1(1,:);% The first intersection.
endnode=pts1(end,:);% The last intersection.
u=sort([firstnode(1),endnode(1)]);% U coordinates of the 2 ending nodes after sorting them.
v=sort([firstnode(2),endnode(2)]);% V coordinates of the 2 ending nodes after sorting them.

% Set the boundary nodes in 4 edges of the parameter domain in
% anticlockwise order.And get their corresponding indices of nodes.
for i=1:2
    indexnodesu{i}=find(nodes(:,1)==i-1);
    nodesu{i}=nodes(indexnodesu{i},:);
    indexnodesv{i}=find(nodes(:,2)==i-1);
    nodesv{i}=nodes(indexnodesv{i},:);
end
[~,temu{1}]=sort(nodesu{1}(:,2),'descend');
[~,temu{2}]=sort(nodesu{2}(:,2));
[~,temv{1}]=sort(nodesv{1}(:,1));
[~,temv{2}]=sort(nodesv{2}(:,1),'descend');
for i=1:2
    indexnodesu{i}=indexnodesu{i}(temu{i});
    nodesu{i}=nodesu{i}(temu{i},:);
    indexnodesv{i}=indexnodesv{i}(temv{i});
    nodesv{i}=nodesv{i}(temv{i},:);
end
nodesuv=[nodesv{1};nodesu{2};nodesv{2};nodesu{1}];
nodesuv=unique(nodesuv,'rows','stable');   
indexnodesuv=[indexnodesv{1};indexnodesu{2};indexnodesv{2};indexnodesu{1}];
indexnodesuv=unique(indexnodesuv,'stable'); 

% The original surface is trimmed into 2 surfaces without holl and there
% are 10 situations of the domain.The 2 boundaries of the 2 trimming
% surfaces are in anticlockwise order.
uv=[u,v];
m=sum(logical((abs(uv-1)<eps | uv<eps)));
if (m==2)
    if (uv(1)<eps)% u1=0
        if (uv(2)<eps)% u1=0,u2=0
            bndpts{1}=nodesu{1}(nodesu{1}(:,2)<uv(4) & nodesu{1}(:,2)>uv(3),:);
            indexbndpts{1}=indexnodesu{1}(nodesu{1}(:,2)<uv(4) & nodesu{1}(:,2)>uv(3),:);
            bndpts{2}=[nodesu{1}((nodesu{1}(:,2)<uv(3)),:);nodesv{1};nodesu{2};nodesv{2};nodesu{1}((nodesu{1}(:,2)>uv(4)),:)];
            indexbndpts{2}=[indexnodesu{1}((nodesu{1}(:,2)<uv(3)),:);indexnodesv{1};indexnodesu{2};indexnodesv{2};indexnodesu{1}((nodesu{1}(:,2)>v(2)),:)];
        elseif (abs(uv(2)-1)<eps)% u1=0,u2=1
            if (u(1)==firstnode(1))
                leftnode=firstnode(2);
                rightnode=endnode(2);
            else
                leftnode=endnode(2);
                rightnode=firstnode(2);
            end
            bndpts{1}=[nodesu{1}(nodesu{1}(:,2)<leftnode,:); nodesv{1};nodesu{2}(nodesu{2}(:,2)<rightnode,:)];
            indexbndpts{1}=[indexnodesu{1}(nodesu{1}(:,2)<leftnode,:); indexnodesv{1};indexnodesu{2}(nodesu{2}(:,2)<rightnode,:)];
            bndpts{2}=[nodesu{2}(nodesu{2}(:,2)>rightnode,:); nodesv{2};nodesu{1}(nodesu{1}(:,2)>leftnode,:)];
            indexbndpts{2}=[indexnodesu{2}(nodesu{2}(:,2)>rightnode,:); indexnodesv{2};indexnodesu{1}(nodesu{1}(:,2)>leftnode,:)];           
        elseif (uv(3)<eps)% u1=0,v1=0          
            bndpts{1}=[nodesu{1}(nodesu{1}(:,2)<uv(4),:);nodesv{1}(nodesv{1}(:,1)<uv(2),:)];
            indexbndpts{1}=[indexnodesu{1}(nodesu{1}(:,2)<uv(4),:);indexnodesv{1}(nodesv{1}(:,1)<uv(2),:)];
            bndpts{2}=[nodesv{1}(nodesv{1}(:,1)>uv(2),:);nodesu{2};nodesv{2};nodesu{1}(nodesu{1}(:,2)>uv(4),:)];
            indexbndpts{2}=[indexnodesv{1}(nodesv{1}(:,1)>uv(2),:);indexnodesu{2};indexnodesv{2};indexnodesu{1}(nodesu{1}(:,2)>uv(4),:)];
        else %u1=0,v2=1
            bndpts{1}=[nodesv{2}(nodesv{2}(:,1)<uv(2),:);nodesu{1}(nodesu{1}(:,2)>uv(3),:)];
            indexbndpts{1}=[indexnodesv{2}(nodesv{2}(:,1)<uv(2),:);indexnodesu{1}(nodesu{1}(:,2)>uv(3),:)];
            bndpts{2}=[nodesu{1}(nodesu{1}(:,2)<uv(3),:);nodesv{1};nodesu{2};nodesv{2}(nodesv{2}(:,1)>uv(2),:)];
            indexbndpts{2}=[indexnodesu{1}(nodesu{1}(:,2)<uv(3),:);indexnodesv{1};indexnodesu{2};indexnodesv{2}(nodesv{2}(:,1)>uv(2),:)];
        end
    elseif (abs(uv(1)-1)<eps)% u1=1,u2=1
        bndpts{1}=nodesu{2}(nodesu{2}(:,2)<uv(4) & nodesu{2}(:,2)>uv(3),:);
        indexbndpts{1}=indexnodesu{2}(nodesu{2}(:,2)<uv(4) & nodesu{2}(:,2)>uv(3),:);
        bndpts{2}=[nodesu{2}((nodesu{2}(:,2)>uv(4)),:);nodesv{2};nodesu{1};nodesv{1};nodesu{2}(nodesu{2}(:,2)<uv(3),:)];
        indexbndpts{2}=[indexnodesu{2}((nodesu{2}(:,2)>uv(4)),:);indexnodesv{2};indexnodesu{1};indexnodesv{1};indexnodesu{2}(nodesu{2}(:,2)<uv(3),:)];
    else % u1~=0,u1~=1
        if (abs(uv(2)-1)<eps)% u2=1
            if (uv(3)<eps)% u2=1,v1=0
                bndpts{1}=[nodesv{1}(nodesv{1}(:,1)>uv(1),:);nodesu{2}(nodesu{2}(:,2)<uv(4),:)];
                indexbndpts{1}=[indexnodesv{1}(nodesv{1}(:,1)>uv(1),:);indexnodesu{2}(nodesu{2}(:,2)<uv(4),:)];                
                bndpts{2}=[nodesu{2}(nodesu{2}(:,2)>uv(4),:);nodesv{2};nodesu{1};nodesv{1}(nodesv{1}(:,1)<uv(1),:)];
                indexbndpts{2}=[indexnodesu{2}(nodesu{2}(:,2)>uv(4),:);indexnodesv{2};indexnodesu{1};indexnodesv{1}(nodesv{1}(:,1)<uv(1),:)];                      
            else % u2=1,v2=1
                bndpts{1}=[nodesu{2}(nodesu{2}(:,2)>uv(3),:);nodesv{2}(nodesv{2}(:,1)>uv(1),:)];
                indexbndpts{1}=[indexnodesu{2}(nodesu{2}(:,2)>uv(3),:);indexnodesv{2}(nodesv{2}(:,1)>uv(1),:)];
                bndpts{2}=[nodesv{2}(nodesv{2}(:,1)<uv(1),:);nodesu{1};nodesv{1};nodesu{2}(nodesu{2}(:,2)<uv(3),:)];
                indexbndpts{2}=[indexnodesv{2}(nodesv{2}(:,1)<uv(1),:);indexnodesu{1};indexnodesv{1};indexnodesu{2}(nodesu{2}(:,2)<uv(3),:)];
            end
        else % u1~=0,u1~=1,u2~=1
            if (uv(3)<eps)% v1=0
                if (uv(4)<eps)% v1=0,v2=0
                    bndpts{1}=nodesv{1}((nodesv{1}(:,1)<uv(2) & nodesv{1}(:,1)>uv(1)),:);
                    indexbndpts{1}=indexnodesv{1}((nodesv{1}(:,1)<uv(2) & nodesv{1}(:,1)>uv(1)),:);
                    bndpts{2}=[nodesv{1}((nodesv{1}(:,1)>uv(2)),:);nodesu{2};nodesv{2};nodesu{1};nodesv{1}((nodesv{1}(:,1)<uv(1)),:)];
                    indexbndpts{2}=[indexnodesv{1}((nodesv{1}(:,1)>uv(2)),:);indexnodesu{2};indexnodesv{2};indexnodesu{1};indexnodesv{1}((nodesv{1}(:,1)<uv(1)),:)];
                else % v1=0,v2=1
                    if (v(1)==firstnode(2))
                        lownode=firstnode(1);
                        upnode=endnode(1);
                    else
                        lownode=endnode(1);
                        upnode=firstnode(1);
                    end
                    bndpts{1}=[nodesv{2}(nodesv{2}(:,1)<upnode,:);nodesu{1};nodesv{1}(nodesv{1}(:,1)<lownode,:)];
                    indexbndpts{1}=[indexnodesv{2}(nodesv{2}(:,1)<upnode,:);indexnodesu{1};indexnodesv{1}(nodesv{1}(:,1)<lownode,:)];
                    bndpts{2}=[nodesv{1}(nodesv{1}(:,1)>lownode,:);nodesu{2};nodesv{2}(nodesv{2}(:,1)>upnode,:)];
                    indexbndpts{2}=[indexnodesv{1}(nodesv{1}(:,1)>lownode,:);indexnodesu{2};indexnodesv{2}(nodesv{2}(:,1)>upnode,:)];
                end
            else % v1=1,v2=1
                bndpts{1}=nodesv{2}((nodesv{2}(:,1)<uv(2) & nodesv{2}(:,1)>uv(1)),:);
                indexbndpts{1}=indexnodesv{2}((nodesv{2}(:,1)<uv(2) & nodesv{2}(:,1)>uv(1)),:);
                bndpts{2}=[nodesv{2}((nodesv{2}(:,1)<uv(1)),:);nodesu{1};nodesv{1};nodesu{2};nodesv{2}((nodesv{2}(:,1)>uv(2)),:)];
                indexbndpts{2}=[indexnodesv{2}((nodesv{2}(:,1)<uv(1)),:);indexnodesu{1};indexnodesv{1};indexnodesu{2};indexnodesv{2}((nodesv{2}(:,1)>uv(2)),:)];
            end
        end
    end
    
    for i=1:2
        bndpts{i}=unique(bndpts{i},'rows','stable');
        indexbndpts{i}=unique(indexbndpts{i},'stable'); %indexbndpts{2}=unique(indexbndpts{i},'stable');
    end
    for iii=1:2
        temp1=norm(bndpts{iii}(1,:)-firstnode);
        temp2=norm(bndpts{iii}(1,:)-endnode);
        if (temp1>temp2)
            boundpts{iii}=[pts1;bndpts{iii};firstnode];
            boundnodesid{iii}=[indexpts1(:);indexbndpts{iii}(:);indexpts1(1)];
        else
            boundpts{iii}=[pts1(end:-1:1,:);bndpts{iii};endnode];
            tem0=indexpts1(end:-1:1);
            tem0=tem0(:);
            boundnodesid{iii}=[tem0;indexbndpts{iii}(:);indexpts1(end)];
        end
    end      
    [in,on]=inpolygon(nodes(:,1),nodes(:,2),boundpts{1}(:,1),boundpts{1}(:,2));
    % MATLAB in=inpolygon, in stands the points which are in the domain or
    %on the boundary of the domain.
    pts{1}=nodes(in,:);
    pts{2}=[pts1;nodes(~in,:)];
%     pts{2}=[boundpts{2}(1:length(pts1),:);nodes(~in,:)];
       
elseif (m==1)
    % The original surface is not trimmed into 2 surface, there's just an
    % intersection curve in the surface and one endding point of the intersection is located in the par-edge.     
    boundpts{1}=[nodesuv;pts1];
    boundnodesid{1}=[indexnodesuv(:);indexpts1(:)];
    pts{1}=nodes;
    in=true(length(nodes),1);
    on=false(length(nodes),1);
    on(boundnodesid{1})=1;
    
else
    %There's a holl in the surface.The outer boundary is in anticlockwise
    %order and the inner boundary is in clockwise order.
    if (norm(firstnode-endnode)<tol)
% Inner loop,the fisrt point is the same as the last one.Determine whether the loop is in clockwise.This can be handled by using
% Green formula.Here the method is based on the cross product in the neighborhood of convex point.
        indexmaxx=find(pts1(:,1)==max(pts1(:,1)));
        if (length(indexmaxx)==1)
            panding1=pts1(indexmaxx-1,:)-pts1(indexmaxx,:);panding1(3)=0;
            panding2=pts1(indexmaxx+1,:)-pts1(indexmaxx,:);panding2(3)=0;
        else
            panding1=pts1(end-1,:)-pts1(end,:);panding1(3)=0;
            panding2=pts(2,:)-pts1(end,:);panding2(3)=0;
        end           
        panding3=cross(panding1,panding2);
        % The cross product is positive when clockwise and negative
        % when anticlockwise.
        if (sum(panding3)<0)
            pts1=pts1(end:-1:1,:);
            indexpts1=indexpts1(end:-1:1,:);
        end               
        boundpts{1}=pts1; 
        boundnodesid{1}=indexpts1(:);
        % outer surface boundary, which consists of both inner and outer loop, notice there is [nan,nan] in boundpts{2} of the boundary of srf2
        boundpts{2}=[nodesuv;nodesuv(1,:);nan,nan;boundpts{1}];
        % The index corresponding to [NAN,NAN] is 0.
        boundnodesid{2}=[indexnodesuv;indexnodesuv(1);0;boundnodesid{1}];
        [in,on]=inpolygon(nodes(:,1),nodes(:,2),boundpts{1}(:,1),boundpts{1}(:,2));
        % in and on are corresponding to surface 1, which in this case
        % is the inner holl.
        pts{1}=nodes(in | on,:);
        pts{2}=nodes(~in | on,:);           
    else            
        %There's just an intersection curve in the surface and the curve is in the inner of parameter domain
        boundpts{1}=[nodesuv;pts1];
        boundnodesid{1}=[indexnodesuv(:);indexpts1(:)];
        pts{1}=nodes; 
        in=true(length(nodes),1);
        on=false(length(nodes),1);
        on(boundnodesid{1})=1;
    end
end

if (nargout==6)
    n=length(boundnodesid);
    boundptsid=cell(1,n);
    if (n==2)
        in_{1}=in;
        in_{2}=~in;
        in_{2}(indexpts1)=1;
    else
        in_{1}=in;
    end
    for i=1:n
        m=length(boundnodesid{i});
        for j=1:m
            tem=boundnodesid{i}(j);
            if (tem==0)
                boundptsid{i}=[boundptsid{i};0];
            else
                boundptsid{i}=[boundptsid{i};sum(in_{i}(1:tem))];
            end
        end
    end
end



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
% k=2;
% figure;hold on;
% triplot(ntri,npts_(:,1),npts_(:,2));
% plot(pts{k}(:,1),pts{k}(:,2),'ro');
% plot(boundpts{k}(:,1),boundpts{k}(:,2),'k','Linewidth',3);
% axis equal;


