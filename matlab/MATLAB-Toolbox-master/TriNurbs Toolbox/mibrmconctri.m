function tr=mibrmconctri(pts1,pts2,id,pnts1,ite)

% Triangulation of 2-D based on the MIBRM(Midpoint insertion boundary recovery
% method) and Delaunay Triangulation. This can make triangulation of an
% arbitrary polygon including concave ones, with all triangular patchs
% in the domain enveloped by boundaries. Notice that this function will
% lead to extra points on the boundaries and NOT delete the triangulars out
% of the boundaries. See conctri.不再使用该函数

% Input:
%   pts1: Parameter coordinates of all the points needed to be triangulated
%   pts2: Parameter coordinates of the points in boundaries. The points have
%       to be sorted in a specific sequence which can be saved in id, usually anticlockwise. The first
%       element is the same as the last one.
%   id: Indices of pts2 in pts1
%   pnts1: Physical coordinates of all the points needed to be triangulated, corresponding to the parameter points pts1
% Output:
%   tr: Triangular structure, including the triangulation connectivity list, the coordinates
%       of points in parameter and physical domain and indices of points on boundaries.

if (nargin==4)   
    ite=0;
end

%修改为含边界限制的剖分，带第3个输入参数，类似delaunaytri
tem1=find(id==0,1);
nn=length(id);
if (isempty(tem1))
    bndid=[id(1:nn-1);id(2:nn-1),id(1)]';
    ppts2=pts2;
    idid=id;
    tem1=0;
else
    tem2=pts2(1,:);
    if (tem2(1)~=0 && tem2(2)~=0 && tem2(1)~=1 && tem2(2)~=1)
        bndid=[id(1:tem1-2);id(2:tem1-2),id(1)]';
        ppts2=pts2(1:tem1-1,:);
        idid=id(1:tem1-1);
    else
        bndid=[id(tem1+1:nn-1);id(tem1+2:nn-1),id(1+tem1)]';
        ppts2=pts2(tem1+1:end,:);
        idid=id(tem1+1:end);
    end
end


TR=delaunayTriangulation(pts1(:,1),pts1(:,2),bndid);
tri=TR.Points;
tricon=TR.ConnectivityList;
%如果遇到自动插入点的问题，也需要调用delaunaydeal函数以解决
tr.delaunay=tricon;
tr.nodes=tri;
% if there is a loop, the index in tr.boundary corresponding to [nan,nan] is 0
tr.boundary=id;
tr.points=pnts1;

n1=length(pts1);%numbers of all the points
n2=length(ppts2);%numbers of points in boundaries,of there is an inner loop, pts2 just restores points of the inner boundary

%,relations between all the points to triangulars, see tnrbpts2tri
nt=length(tricon);%numbers of triangulars
p2t=cell(1, n1);
for ii=1:nt
    for jj=1:3
        p2t{tricon(ii,jj)}=[p2t{tricon(ii,jj)}, ii];
    end
end

for i=1:n2-1
    %relations between points to edges, see tnrbpt2edges
    pt=idid(i);
    temtri=p2t{pt};
    n=length(temtri);
    edges=zeros(2*n, 1);
    for ii=1:n
        p=tricon(temtri(ii),:)~=pt;
        edges(2*ii-1:2*ii)=tricon(temtri(ii),p);
    end
    edges=RemDuplicate(edges);
    
    %whether the boundary edge is in the union of triangulation edges
    de=logical(edges==idid(i+1));
    if (sum(de)==0)
        if (ite==5) %MIBRM无法去除插入点，改用其他算法；最多插入5次就强制停止
            continue;
        else
            cpoint=(ppts2(i,:)+ppts2(i+1,:))/2;
            pts2(tem1+i+2:end+1,:)=pts2(tem1+i+1:end,:);
            pts2(tem1+i+1,:)=cpoint;
            pts1(end+1,:)=cpoint;
            cpoint_=(pnts1(idid(i),:)+pnts1(idid(i+1),:))/2;
            pnts1(end+1,:)=cpoint_;
            
            id(tem1+i+2:end+1)=id(tem1+i+1:end);
            id(tem1+i+1)=length(pts1);
            ite=ite+1;
            
            %Recursive call the function
            tr=mibrmconctri(pts1,pts2,id,pnts1,ite);
        end
        break;
    end
end
        
        


