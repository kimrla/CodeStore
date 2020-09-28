function [sed1, stri2, spts1, spts2, spnts1, spnts2]=tnrbintersects(tnrb1, tnrb2, p2t1, p2t2)

% tnrbintersects: Get the intersection points of two tri-nurbs surfaces and sort them
% 
% Calling Sequences:
% 
%       [sed1, stri2, spts1, spts2, spnts1, spnts2]=tnrbintersects(tnrb1, tnrb2, p2t1, p2t2)
% 
% INPUTS:
% 
%       tnrb1, tnrb2 - Triangular representation of two nurbs surface.
% 
%       p2t1, p2t2 - The relations from points to triangles of two 
%                     tri-nurbs surface. See also tnrbpts2tri.
%
% OUTPUT: 
%
%      sed1  -  Edges of the triangles of tri-nurbs surface 1 that
%                   intersected with tri-nurbs surface 2.
%
%      stri2   -  Triangles tri-nurbs surface 2 that intersected with 
%                    tri-nurbs surface 1.
%
%       spts1, spts2 - Parametric intersection points of the two surfaces.
%
%       spnts1, spnts2 - Intersection points of the two surfaces.
%
%  Discription:
%   
%      This routine used the 3D triangles to get interections.
%

% The nearest points of the two surfaces
tol=max([tnrb1.seeds(1), tnrb2.seeds(1)]);
[p1, p2]=nearpnts(tnrb1.points, tnrb2.points, tol);

% Intersections of two tri-nurbs surfaces
[ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri(tnrb1, tnrb2, p2t1, p2t2, p1, p2);

% Sort the intersections of two surfaces and solve losted interections
r=1; ne=size(ed1,1);
while ne>1 && r<100
    % Get a new sorted interections
    [sq, id]=tnrbintersort(tnrb1, ed1, pts1, p2t1);%经过该函数排好序的点只是所有交点的一部分，可能有些点未参与排序
    %通过不断循环，将所有交点排好序，只不过拍完后交点可能被分割成多个分支（r条），所有交点没有连成一条曲线
    % Add to sorted data
    sed1{r}=ed1(sq,:);
    stri2{r}=tri2(sq,:);
    spts1{r}=pts1(sq,:);
    spts2{r}=pts2(sq,:);
    spnts1{r}=pnts1(sq,:);
    spnts2{r}=pnts2(sq,:);    
    
    % The remained intersections
    ed1=ed1(id,:);
    tri2=tri2(id,:);
    pts1=pts1(id,:);
    pts2=pts2(id,:);
    pnts1=pnts1(id,:);
    pnts2=pnts2(id,:);
    ne=size(ed1,1);  
    
%     ne=sum(id);
    
    [sed1{r}, id]=RemDuplicate(sed1{r});
    stri2{r}=stri2{r}(id,:);
    spts1{r}=spts1{r}(id,:);
    spts2{r}=spts2{r}(id,:);
    spnts1{r}=spnts1{r}(id,:);
    spnts2{r}=spnts2{r}(id,:);
    r=r+1;
end

if r==1
    sed1{r}=[];
    stri2{r}=[];
    spts1{r}=[];
    spts2{r}=[];
    spnts1{r}=[];
    spnts2{r}=[];
end

% Connect and extend the intersections
nr=length(sed1); r=1;
while nr>=1 && r<=nr
    % Search from beginning
    do=true; t=1;
    while do && t<100
        ei1=sed1{r}(1,:);
        trg0=sort(RemDuplicate([sed1{r}(1,:), sed1{r}(2,:)]'))';%trg0对应的三角形其中两条边含有交点
        tri=tnrbedge2tri(p2t1, ei1); 
        tri=tnrb1.delaunay(tri,:);
        d=DistanceMatrix(trg0, tri);
        trg1=tri(d~=0,:);%trg1对应的三角形其中1条边含有交点，另外两条边是利用tnrbintertri2tri待确定的
        trg2=stri2{r}(1,:);
        if isempty(trg1)
            break;
        else
            [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri2tri(tnrb1, tnrb2, p2t2, ei1, trg1, trg2);%调用tnrbintertri2tri子函数
            if isempty(ed1)
                break;
            end
            ad=false;
            s=1;
            while nr>=1 && s<=nr
                d=DistanceMatrix(ed1, sed1{s});
                nd=find(d==0);
                if isempty(nd)
                    nd=0;
                elseif length(nd)>1
                    break;
                end
                ns=length(sed1{s});
                if nd==ns
                    if s==r
                        sed1{r}=[ed1; sed1{r}]; 
                        %如果通过tnrbintertri2tri能找到新的交点，则将其作为第一个点放在已得到的交线分支上（因为是从beginning点开始查找的）
                        %并更新所有项，将新求得的点作为该分支上首个交点继续进行循环
                        stri2{r}=[tri2; stri2{r}]; 
                        spts1{r}=[pts1; spts1{r}]; 
                        spts2{r}=[pts2; spts2{r}]; 
                        spnts1{r}=[pnts1; spnts1{r}]; 
                        spnts2{r}=[pnts2; spnts2{r}]; 
                        ad=true;
                        break;
                    else
                        sed1{r}=[sed1{s}; sed1{r}]; 
                        stri2{r}=[stri2{s}; stri2{r}]; 
                        spts1{r}=[spts1{s}; spts1{r}]; 
                        spts2{r}=[spts2{s}; spts2{r}]; 
                        spnts1{r}=[spnts1{s}; spnts1{r}]; 
                        spnts2{r}=[spnts2{s}; spnts2{r}]; 
                        sed1(s)=[]; 
                        stri2(s)=[]; 
                        spts1(s)=[]; 
                        spts2(s)=[]; 
                        spnts1(s)=[]; 
                        spnts2(s)=[]; 
                        nr=nr-1;
                        ad=true;
                    end                    
                elseif nd==1 && s~=r
                    p=ns:-1:1;
                    sed1{r}=[sed1{s}(p,:); sed1{r}]; 
                    stri2{r}=[stri2{s}(p,:); stri2{r}]; 
                    spts1{r}=[spts1{s}(p,:); spts1{r}]; 
                    spts2{r}=[spts2{s}(p,:); spts2{r}]; 
                    spnts1{r}=[spnts1{s}(p,:); spnts1{r}]; 
                    spnts2{r}=[spnts2{s}(p,:); spnts2{r}]; 
                    sed1(s)=[]; 
                    stri2(s)=[]; 
                    spts1(s)=[]; 
                    spts2(s)=[]; 
                    spnts1(s)=[]; 
                    spnts2(s)=[]; 
                    nr=nr-1;
                    ad=true;
                end
                s=s+1;
            end
            if ad==false
                sed1{r}=[ed1; sed1{r}]; 
                stri2{r}=[tri2; stri2{r}]; 
                spts1{r}=[pts1; spts1{r}]; 
                spts2{r}=[pts2; spts2{r}]; 
                spnts1{r}=[pnts1; spnts1{r}]; 
                spnts2{r}=[pnts2; spnts2{r}]; 
            end
        end
        t=t+1;
    end
    
    % Search from end
    do=true; t=1;
    while do && t<100
        ei1=sed1{r}(end,:);
        trg0=sort(RemDuplicate([sed1{r}(end-1,:), sed1{r}(end,:)]'))';
        tri=tnrbedge2tri(p2t1, ei1); 
        tri=tnrb1.delaunay(tri,:);
        d=DistanceMatrix(trg0, tri);
        trg1=tri(d~=0,:);
        trg2=stri2{r}(end,:);
        if isempty(trg1)
            break;
        else
            [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri2tri(tnrb1, tnrb2, p2t2, ei1, trg1, trg2);
            if isempty(ed1)
                break;
            end
            ad=false;
            s=1;
            while nr>=1 && s<=nr
                d=DistanceMatrix(ed1, sed1{s});
                nd=find(d==0);
                if isempty(nd)
                    nd=0;
                elseif length(nd)>1
                    break;
                end
                ns=length(sed1{s});
                if nd==1
                    if s==r
                        sed1{r}=[sed1{r}; ed1]; 
                        stri2{r}=[stri2{r}; tri2]; 
                        spts1{r}=[spts1{r}; pts1]; 
                        spts2{r}=[spts2{r}; pts2]; 
                        spnts1{r}=[spnts1{r}; pnts1]; 
                        spnts2{r}=[spnts2{r}; pnts2]; 
                        ad=true;
                        break;
                    else
                        sed1{r}=[sed1{r}; sed1{s}]; 
                        stri2{r}=[stri2{r}; stri2{s}]; 
                        spts1{r}=[spts1{r}; spts1{s}]; 
                        spts2{r}=[spts2{r}; spts2{s}]; 
                        spnts1{r}=[spnts1{r}; spnts1{s}]; 
                        spnts2{r}=[spnts2{r}; spnts2{s}]; 
                        sed1(s)=[]; 
                        stri2(s)=[]; 
                        spts1(s)=[]; 
                        spts2(s)=[]; 
                        spnts1(s)=[]; 
                        spnts2(s)=[]; 
                        nr=nr-1;
                        ad=true;
                    end                    
                elseif nd==ns && r~=s
                    p=ns:-1:1;
                    sed1{r}=[sed1{r}; sed1{s}(p,:)]; 
                    stri2{r}=[stri2{r}; stri2{s}(p,:)]; 
                    spts1{r}=[spts1{r}; spts1{s}(p,:)]; 
                    spts2{r}=[spts2{r}; spts2{s}(p,:)]; 
                    spnts1{r}=[spnts1{r}; spnts1{s}(p,:)]; 
                    spnts2{r}=[spnts2{r}; spnts2{s}(p,:)]; 
                    sed1(s)=[]; 
                    stri2(s)=[]; 
                    spts1(s)=[]; 
                    spts2(s)=[]; 
                    spnts1(s)=[]; 
                    spnts2(s)=[]; 
                    nr=nr-1;
                    ad=true;
                end
                s=s+1;
            end
            if ad==false
                sed1{r}=[sed1{r}; ed1]; 
                stri2{r}=[stri2{r}; tri2]; 
                spts1{r}=[spts1{r}; pts1]; 
                spts2{r}=[spts2{r}; pts2]; 
                spnts1{r}=[spnts1{r}; pnts1]; 
                spnts2{r}=[spnts2{r}; pnts2]; 
            end
        end
        t=t+1;
    end
    r=r+1;
end

%% demo
% % The mesh seed length (h0)
% h0=0.5;
% 
% % Create a nurbs sphere
% center=[8,5,2];
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
% 
% figure; hold on;
% triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2));  
% for r=1:nr
%     plot(spts1{r}(:,1), spts1{r}(:,2), 'k.', 'MarkerSize', 13); 
%     plot(spts1{r}(:,1), spts1{r}(:,2), 'r', 'LineWidth', 1); 
% end
% title('Parametric mesh of surface 1');  
% axis equal; 
% 
% figure; hold on; 
% triplot(tnrb2.delaunay, tnrb2.nodes(:,1), tnrb2.nodes(:,2)); 
% for r=1:nr
%     plot(spts2{r}(:,1), spts2{r}(:,2), 'k.', 'MarkerSize', 13); 
%     plot(spts2{r}(:,1), spts2{r}(:,2), 'r', 'LineWidth', 1); 
% end
% title('Parametric mesh of surface 2'); 
% axis equal;







