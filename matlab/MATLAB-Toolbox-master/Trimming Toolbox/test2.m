TOL=5*10e-3;
quad=[];
TRI=[];

tri=newtri;
pts_now=pts{k};
pnts_now=nnpnts{k};

e1=tri(:, [1,2]);
e2=tri(:, [2,3]);
e3=tri(:, [3,1]);
ed=sort([e1; e2; e3], 2);
ed=unique(ed,'rows');

p2t=cell(1,length(pts_now));%细胞数组
for i=1:length(tri)
    for j=1:3
        p2t{tri(i,j)}=[p2t{tri(i,j)}, i];
    end
end

for i=1:length(ed)
    tri1=p2t{ed(i,1)};
    tri2=p2t{ed(i,2)};
    tri3=[tri1, tri2];
    [~, id]=RemDuplicate(tri3');
    tri3=tri3(~id);%矢量，非细胞数组
    if length(tri3)==2
        p3=setdiff(tri(tri3(1),:),ed(i,:));
        p4=setdiff(tri(tri3(2),:),ed(i,:));
        if abs(norm(pts_now(ed(i,1),:)-pts_now(ed(i,2),:))-norm(pts_now(p3,:)-pts_now(p4,:)))<TOL
            quad=[quad;[ed(i,1),p3,ed(i,2),p4]];         
            TRI=[TRI,tri3];
        end       
    end
  
end

tri(TRI,:)=[];
tri(:,4)=NaN;
newtri_now=[tri;quad];

% figure;hold on;
% quadmesh(newtri,pts_now(:,1),pts_now(:,2));
% % plot(boundpts{k}(:,1),boundpts{k}(:,2),'r--','Marker','*','MarkerSize',4);
% axis equal;

figure;hold on;
quadsurf(newtri_now,pnts_now(:,1),pnts_now(:,2),pnts_now(:,3));
view(3);
axis equal;