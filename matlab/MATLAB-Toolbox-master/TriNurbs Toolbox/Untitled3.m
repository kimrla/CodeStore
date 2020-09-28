h0=0.8;

% Create a nurbs sphere
center=[5,5,4];
circ=nrbcirc(4, center, 0, pi);
srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
srf2=nrbtestsrf;

% Transform a nurbs surface into triangular representation
tnrb1=nrb2tri(srf1, h0); 
tnrb2=nrb2tri(srf2, h0); 

% Get the edges of a tri-nurbs surface that intersected with another tri-nurbs surface
p2t1=tnrbpts2tri(tnrb1);
p2t2=tnrbpts2tri(tnrb2);

% Get the intersection points of two tri-nurbs surfaces and sort them
[sed1, stri2, spts1, spts2, spnts1, spnts2]=tnrbintersects(tnrb1, tnrb2, p2t1, p2t2);

% Connect and extend the intersections
inter1=tnrbinterconnct(tnrb1, tnrb2, p2t1, p2t2, {sed1, stri2, spts1, spts2, spnts1, spnts2});
sed1=inter1{1}; stri2=inter1{2};
spts1=inter1{3}; spts2=inter1{4};
spnts1=inter1{5}; spnts2=inter1{6};

cc=find(spts1{1}(:,1)==spts1{1}(2,1) & spts1{1}(:,2)==spts1{1}(2,2));
if (length(cc)~=1)
    cc=cc(2);
    sed1{1}(cc:end,:)=[];
    stri2{1}(cc:end,:)=[];
    spts1{1}(cc:end,:)=[];
    spts2{1}(cc:end,:)=[];
    spnts1{1}(cc:end,:)=[];
    spnts2{1}(cc:end,:)=[];
end

dist=h0;
[npts,npnts1,npnts2,ndist]=equinter(spnts1{1},tnrb1,dist);
tritol=h0/4;
[npts_,npnts_,ntri,intid]=ndelaunaytri(tnrb1,npts,npnts2,tritol);

indexpts1=[intid(:,1);intid(end,2)]';
[pts,in,on,boundpts,boundnodesid,boundptsid]=npointinparsrf(npts_,indexpts1);

nnpnts=cell(1,2);
nnpnts{1}=npnts_(in,:);
% nnpnts{2}=[npnts_(indexpts1,:);npnts_(~in,:)]; % 没有挖孔，把交点放在最前面
nnpnts{2}=npnts_(~in | on,:); % 有挖孔

% k=2; 
% in_=~in;
% in_(indexpts1)=1; 

k=1;
in_=in;

newtri=transindex(ntri,in_); % newtri可能在有效域外还有三角形
newtri_=nconctri(newtri,pts{k},boundptsid{k});

% figure;hold on;
% triplot(newtri_.delaunay,pts{k}(:,1),pts{k}(:,2));
% % plot(boundpts{k}(:,1),boundpts{k}(:,2),'r--','Marker','*','MarkerSize',4);
% axis equal;

figure;hold on;
triplot(newtri,pts{k}(:,1),pts{k}(:,2));
% plot(boundpts{k}(:,1),boundpts{k}(:,2),'r--','Marker','*','MarkerSize',4);
axis equal;

% figure;hold on;
% trisurf(newtri_.delaunay,nnpnts{k}(:,1),nnpnts{k}(:,2),nnpnts{k}(:,3));
% axis equal;

figure;hold on;
trisurf(newtri,nnpnts{k}(:,1),nnpnts{k}(:,2),nnpnts{k}(:,3));
axis equal;