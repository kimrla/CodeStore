function inter=tnrbinterconnct(tnrb1, tnrb2, p2t1, p2t2, inter)

% tnrbinterconnct: Connect and extend the intersections
% 
% Calling Sequences:
% 
%       inter=tnrbinterconnct(tnrb1, tnrb2, p2t1, p2t2, inter)
% 
% INPUTS:
% 
%       tnrb1, tnrb2 - Triangular representation of two nurbs surface.
% 
%       p2t1, p2t2 - The relations from points to triangles of two 
%                     tri-nurbs surface. See also tnrbpts2tri.
%
%       inter = {sed1, stri2, spts1, spts2, spnts1, spnts2}
%                   The outputs of tnrbintersects.
%
% OUTPUT: 
%
%      inter = {sed1, stri2, spts1, spts2, spnts1, spnts2}
%                  The intersections after connection and extension.
%

sed1=inter{1}; stri2=inter{2};
spts1=inter{3}; spts2=inter{4};
spnts1=inter{5}; spnts2=inter{6};

nr=length(sed1); r=1;
while nr>=1 && r<=nr
    % Search from beginning
    do=true; t=1;
    while do && t<100
        ei1=sed1{r}(1,:);
        trg0=sort(RemDuplicate([sed1{r}(1,:), sed1{r}(2,:)]'))';
        tri=tnrbedge2tri(p2t1, ei1); 
        tri=tnrb1.delaunay(tri,:);
        d=DistanceMatrix(trg0, tri);
        trg1=tri(d~=0,:);
        trg2=stri2{r}(1,:);
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
                if nd==ns
                    if s==r
                        sed1{r}=[ed1; sed1{r}]; 
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
inter={sed1, stri2, spts1, spts2, spnts1, spnts2};


%% demo
% % The mesh seed length (h0)
% h0=0.8;
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
% inter=tnrbinterconnct(tnrb1, tnrb2, p2t1, p2t2, {sed1, stri2, spts1, spts2, spnts1, spnts2});
% sed1=inter{1}; stri2=inter{2};
% spts1=inter{3}; spts2=inter{4};
% spnts1=inter{5}; spnts2=inter{6};
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





