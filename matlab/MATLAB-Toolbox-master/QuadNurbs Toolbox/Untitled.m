%demo
n=10;
x=1:n;
y=1:n;
[X,Y]=meshgrid(x,y);

bnd=[12,33,32,53,54,57,15,64];
bnd=[bnd(1:end);bnd(2:end),12]';

tr=delaunayTriangulation(X(:),Y(:),bnd);

pnt=tr.Points;
tri=tr.ConnectivityList;
cnt=tr.Constraints;


m=length(cnt);
figure;
hold on;
triplot(tri,pnt(:,1),pnt(:,2));
plot(pnt(bnd,1),pnt(bnd,2),'ro','Markersize',8);
for i=1:m
    plot(pnt(cnt(i,:),1),pnt(cnt(i,:),2),'dg','LInewidth',2,'Markersize',5);
end

nbndin=delaunaydeal(pnt,cnt,bnd,bnd);