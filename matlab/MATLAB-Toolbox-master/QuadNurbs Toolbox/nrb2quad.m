function qnrb=nrb2quad(nrbs, h0)
%
% nrb2quad: Transform a nurbs surface into quadrangular representation.
% 
% Calling Sequences:
% 
%     qnrb=nrb2quad(srf, h0)
% 
% INPUTS:
%
%      nurbs - A nurbs surface.
%
%      h0 - Mesh seed length.
%
% OUTPUT:
% 
%     qnrb - Quadrangular representation of the nurbs (quad-nurbs) surface.
%               A structure array with the following fields:
%
%         qnrb.form =  'Quad-NURBS'.  Quadrangular representation 
%                               of nurbs surface
%
%         qnrb.surface - The nurbs surface.
%
%         qnrb.numbers - The number of nodes or points
%                                  and the number of quadrangles.
%
%         qnrb.seeds    -  Mesh seed length of the surface.
%
%         qnrb.grids    -  Grids of two directions on parametric domain.
%
%         qnrb.nodes  -  Nodes on parametric domain.
%
%         qnrb.points  -  Points on the nurbs surface.
%
%         qnrb.quad  -  Quadrangular reprezetation.
%
%         qnrb.edges  -  Edges of the quadrangles.
%
%         qnrb.ed2qd  -  Edges to quadrangles.
%
%         qnrb.pt2qd  -  Points to quadrangles.
%
%         qnrb.tri  - Triangulation of the quadrangles. Triangles that
%                          are acturally a line are removed.
%  
%         qnrb.pt2tri - The relations of points with triangles.
%
%         qnrb.UniqPt - Points without dumpilation.
%
%         qnrb.UniqId={nums, dump, uniq}
%
%                nums - Global node numbers of each node. 
%                            qnrb.points=qnrb.UniqPt(nums, :).
%
%                dump -  Indexes correspond to each unique node. 
%                               dump{k} includes all dumplicated nodes at qnrb.UniqPt(k, :).
%
%                uniq  -   Indexes of unique points. qnrb.UniqPt=qnrb.points(uniq, :).
%  

% Prepare nodes and points
[Lx, Ly]=nrbsrfmeasure(nrbs);
m=round(max(Ly)/h0)+1; 
n=round(max(Lx)/h0)+1;
s=linspace(0, 1, m);
t=linspace(0, 1, n);
[y, x]=meshgrid(t, s);
nodes=[x(:), y(:)];
points=nrbeval(nrbs, {s, t});
points=[points(1,:)', points(2,:)', points(3,:)'];

% Quad mesh
quad=zeros((m-1)*(n-1), 4);
for j=1:n-1
    for i=1:m-1
        p=(j-1)*(m-1)+i;
        q=(j-1)*m+i;
        quad(p, 1)=q;
        quad(p, 2)=q+1;
        quad(p, 3)=q+m+1;
        quad(p, 4)=q+m;
    end
end

% Quad edges and their relations with quadrangles
qedges=zeros(m*(n-1)+n*(m-1), 2);
ed2qd=cell(m*(n-1)+n*(m-1), 1);
for j=1:n
    for i=1:m-1
        p=(j-1)*(m-1)+i;
        q=(j-1)*m+i;
        qedges(p,1)=q; 
        qedges(p,2)=q+1; 
        if j==1
            ed2qd{p}=i;
        elseif j==n
            ed2qd{p}=(n-2)*(m-1)+i;
        else
             ed2qd{p}=[(j-2)*(m-1)+i, (j-1)*(m-1)+i];
        end
    end
end
for j=1:n-1
    for i=1:m
        p=(j-1)*m+i;
        qedges(n*(m-1)+p,1)=p; 
        qedges(n*(m-1)+p,2)=p+m; 
        if i==1
            ed2qd{n*(m-1)+p}=(j-1)*(m-1)+i;
        elseif i==m
            ed2qd{n*(m-1)+p}=(j-1)*(m-1)+m-1;
        else
            ed2qd{n*(m-1)+p}=[(j-1)*(m-1)+i-1, (j-1)*(m-1)+i];
        end
    end
end

% Points and their relations with quadrangles
pt2qd=cell(m*n, 1);
for j=1:n
    for i=1:m
        p=(j-1)*m+i;
        if i==1 && j==1
            pt2qd{p}=1;
        elseif i==m && j==1
            pt2qd{p}=m-1;
        elseif i==1 && j==n
            pt2qd{p}=(n-2)*(m-1)+i;
        elseif i==m && j==n
            pt2qd{p}=(n-2)*(m-1)+i-1;
        elseif i>1 && i<m && j==1
            pt2qd{p}=[(j-1)*(m-1)+i-1, (j-1)*(m-1)+i];
        elseif i>1 && i<m && j==n
            pt2qd{p}=[(j-2)*(m-1)+i-1, (j-2)*(m-1)+i];
        elseif i==1 && j>1 && j<n
            pt2qd{p}=[(j-2)*(m-1)+i, (j-1)*(m-1)+i];
        elseif i==m && j>1 && j<n
            pt2qd{p}=[(j-2)*(m-1)+i-1, (j-1)*(m-1)+i-1];
        else
            pt2qd{p}=[(j-2)*(m-1)+i-1, (j-1)*(m-1)+i-1, (j-2)*(m-1)+i, (j-1)*(m-1)+i];
        end
    end
end

% Transform the quad mesh to triangular mesh
tri=[quad(:,[1,2,3]); quad(:,[1,3,4])];

% Remove triangles that are a line
d1=points(tri(:,1),:)-points(tri(:,2),:);
d2=points(tri(:,1),:)-points(tri(:,3),:);
d3=points(tri(:,3),:)-points(tri(:,2),:);
d1=sqrt(sum(d1.^2, 2));
d2=sqrt(sum(d2.^2, 2));
d3=sqrt(sum(d3.^2, 2));
p=d1>h0*1e-6;
p=p & d2>h0*1e-6;
p=p & d3>h0*1e-6;
tri=tri(p,:);

% Points and their relations with triangles
pt2tri=cell(m*n, 1);
q=zeros(size(p));
tt=1;
for i=1:length(p)
    if p(i)
        q(i)=tt; tt=tt+1;
    end
end
for p=1:m*n
    tp=q([pt2qd{p}, (m-1)*(n-1)+pt2qd{p}]);
    id=tp~=0;
    pt2tri{p}=tp(id)';
end

% Find the unique points
[nums, UniqPt, dump, uniq]=numbering(points);

% Make quadrangular representation of nurbs surface
qnrb.form='Quad-NURBS';
qnrb.dim=2;
qnrb.nurbs=nrbs;
qnrb.numbers=[m, n, (m-1)*(n-1)];
qnrb.seeds=h0;
qnrb.grids={s, t};
qnrb.nodes=nodes;
qnrb.points=points;
qnrb.quad=quad;
qnrb.qedges=qedges;
qnrb.ed2qd=ed2qd;
qnrb.pt2qd=pt2qd;
qnrb.tri=tri;
qnrb.pt2tri=pt2tri;
qnrb.UniqPt=UniqPt;
qnrb.UniqId={nums, dump, uniq};


%% demo
% % The mesh seed length (h0)
% h0=1.1;
% 
% % Create a nurbs sphere
% center=[2,2,4];
% circ=nrbcirc(4, center, 0, pi);
% srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% srf=srf2;
% 
% % Transform a nurbs surface into quadrangular representation
% qnrb=nrb2quad(srf, h0);
% points=qnrb.points;
% nodes=qnrb.nodes;
% edges=qnrb.qedges;
% ed2qd=qnrb.ed2qd;
% pt2qd=qnrb.pt2qd;
% 
% % Plot parametric plane
% figure;
% quadplot(qnrb, 'r');
% view(2); axis equal;
% 
% % Plot results
% figure; hold on;
% quadsurf(qnrb.quad, points(:,1), points(:,2), points(:,3));
% view(3); axis equal;
% 
% % Plot an edge and related quadrangles
% k=200;
% plot3(points(edges(k,:),1), points(edges(k,:),2), points(edges(k,:),3), 'LineWidth', 2, 'Color', 'r');
% quadmesh(qnrb.quad(ed2qd{k},:), points(:,1), points(:,2), points(:,3));
% 
% % Plot a points and related quadrangles
% k=102;
% plot3(points(k,1), points(k,2), points(k,3), 'r.', 'MarkerSize', 20);
% quadmesh(qnrb.quad(pt2qd{k},:), points(:,1), points(:,2), points(:,3));






