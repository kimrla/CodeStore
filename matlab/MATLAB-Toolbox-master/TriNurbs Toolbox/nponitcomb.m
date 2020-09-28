function [npts,npnts]=nponitcomb(pts1,pts2,pnts1,pnts2,sed,stri_)
% Combine the 2 series of intersections obtained by tnrbinterconnct or tnrbintersort into 1
% series intersections.

% Input:
%   pts1: Parameter coordinates of intersections on the grids.
%   pts2: Parameter coordinates of intersections which is NOT on the grids
%       but inside the triangulars.
%   pnts1: Coordinates of intersections on the grids.
%   pnts2: Coordinates of intersections which is NOT on the grids
%       but inside the triangulars.
%   sed: Edges of intersections of tnrb1.
%   stri_: Triangulars of intersections of tnrb1.

% Output:
%   npts: New parameter intersections combined by 2 series original points.
%   npnts: New coordinates of intersections combined by 2 series original points.


m=length(pnts1);
n=length(pnts2);
index=1;
tem1=1;
npts=pts1;
npnts=pnts1;
det1=0;
stri=stri_;
intsed=sed;

stri=sort(stri,2);
[mm,~,v]=unique(stri,'rows');
[v,indv]=sort(v);
stri=mm(v,:);

pts2=pts2(indv,:);
pnts2=pnts2(indv,:);

for i=2:n
    tem4=sort(stri(i,:))==sort(stri(i-1,:));
    qunima=sum(tem4);
    if (qunima~=3)
        temtri{tem1}=index:i-1;
        tem1=tem1+1;
        index=i;
    end
end
temtri{tem1}=index:n;
        
numtemtri=length(temtri);
for i=1:numtemtri
    numpts=length(temtri{i});
    for j=1:m-1
        log1=logical(det1==j);
        detsum=sum(log1);
        if (detsum==0)
            cp=unique([intsed(j,:),intsed(j+1,:)]);
            if (cp==sort(stri(temtri{i}(1),:)))
                npts(j+1+numpts:end+numpts,:)=npts(j+1:end,:);
                npnts(j+1+numpts:end+numpts,:)=npnts(j+1:end,:);
                intsed(j+1+numpts:end+numpts,:)=intsed(j+1:end,:);
                intpts=pts2(temtri{i},:);
                intpnts=pnts2(temtri{i},:);
                
                nn=numpts;
                firstpt=npts(j,:);
                insintpts=intpts;
                insintpnts=intpnts;
                newintpts=zeros(numpts,2);
                newintpnts=zeros(numpts,3);
                for p=1:numpts-1
                    distpts=zeros(1,nn);
                    for k=1:nn
                        distpts(k)=norm(firstpt-insintpts(k,:));
                    end
                    [~,pv]=sort(distpts);
                    newintpts(p,:)=insintpts(pv(1),:);
                    newintpnts(p,:)=insintpnts(pv(1),:);
                    nn=nn-1;
                    firstpt=newintpts(p,:);
                    insintpts(pv(1),:)=[];
                    insintpnts(pv(1),:)=[];
                end
                newintpts(end,:)=insintpts(1,:);
                newintpnts(end,:)=insintpnts(1,:);
                intpts=newintpts;
                intpnts=newintpnts;
                               
%                 detx=pts1(j+1,1)-pts1(j,1);
%                 detintx=intpts(end,1)-intpts(1,1);
%                 %按照u向参数排序，只考虑一般情况，特殊情况的交点（如斜率无穷大处附近的交点）可能出现错误
%                 if (detx*detintx<0)
%                     intpts=intpts(end:-1:1,:);
%                     intpnts=intpnts(end:-1:1,:);
%                 elseif( detx*detintx==0)
%                     if (detintx~=0)
%                         dety=pts1(j+1,2)-pts1(j,2);
%                         detinty=intpts(end,2)-intpts(1,2);
%                         if (dety*detinty<0)
%                             intpts=intpts(end:-1:1,:);
%                             intpnts=intpnts(end:-1:1,:);
%                         end
%                     end                   
%                 end

                npts(j+1:j+numpts,:)=intpts;
                npnts(j+1:j+numpts,:)=intpnts;
                intsed(j+1:j+numpts,:)=zeros(numpts,2);
                det1=[det1,j];
            end
        end
    end
    
end
% % The mesh seed length (h0)
% h0=1.8;
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
% [input1,shanchu,~]=unique(spts2_{1},'rows');
% input2=stri2_{1}(shanchu,:);
% input3=spnts2_{1}(shanchu,:);
% [npts,npnts]=nponitcomb(spts1{1},spts2_{1},spnts1{1},spnts2_{1},sed1{1},stri2_{1});
% 
% 
% nds=tnrb1.nodes;
% nr=length(spnts1);
% 
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
% 
% 
% pm=length(npts);
% figure; hold on; 
% triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2)); 
% axis equal;
% plot(npts(:,1), npts(:,2), 'k.', 'MarkerSize', 13); 
% plot(npts(:,1), npts(:,2), 'r', 'LineWidth', 1); 








