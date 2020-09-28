n=10;
x=1:n;
y=1:n;
[X,Y]=meshgrid(x,y);

bnd=[12,33,32,53,54,57,15];
bnd=[bnd(1:end);bnd(2:end),64]';

tr=delaunayTriangulation(X(:),Y(:),bnd);

pnt=tr.Points;
tri=tr.ConnectivityList;
cnt=tr.Constraints;


ncnt=[12,33,32,53,54,101,55,56,57,36,15,101];
ncnt=[ncnt;ncnt(2:end),64]';
tr1=delaunayTriangulation(pnt,ncnt);
pnt1=tr1.Points;
tri1=tr1.ConnectivityList;
cnt1=tr1.Constraints;



m=length(cnt);
figure;
hold on;
triplot(tri,pnt(:,1),pnt(:,2));
plot(pnt(bnd,1),pnt(bnd,2),'ro','Markersize',8);
for i=1:m
    plot(pnt(cnt(i,:),1),pnt(cnt(i,:),2),'dg','LInewidth',2,'Markersize',5);
end


n=length(cnt);
figure;
hold on;
triplot(tri1,pnt1(:,1),pnt1(:,2));
plot(pnt1(ncnt,1),pnt1(ncnt,2),'ro','Markersize',8);
for i=1:n
    plot(pnt1(cnt1(i,:),1),pnt1(cnt1(i,:),2),'dg','LInewidth',2,'Markersize',5);
end
