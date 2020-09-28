function [TRI,PT]=DelaunayWithGrid(Poly,dl,PolyIn)
% function:	DelaunayMeshGrid
% Description:	delaunay triangulation for 2D geometry. geometry can have holes inside 
%               Poly is essential while dl and PolyIn can be vacant
%              
% Input:  
%    Poly   -  n*2,  polygon coordinate of outer layer
%    dl     - approximate element size
%    PolyIn - cell array, each one is the polygon coordinate of inner layer
% Output: 
%    TRI    - Each row of TRI specifies a triangle defined by indices with respect to the points.
%    PT     - m*2, points array
%
%    Example:
%    Poly=[-10 -10;10 -10;10 10;-10 10;];
%    PolyIn={[-7 2;-3 2;-3 6;-7 6],[3 2;7 2;7 6;3 6],[0 -6;3 -3;-3 -3]};
%    [TRI,PT]=DelaunayWithGrid(Poly,0.37,PolyIn);
%
% Routines: p_poly_dist.m
% Revision history:
%    8/6/2013  - by Liufeng  
%                email: liuf_hb@sina.com

if nargin==1
    dl=[];PolyIn=[];
elseif nargin==2
    PolyIn=[];
end

% delete the  repeated points
if Poly(end,1)==Poly(1,1) && Poly(end,2)==Poly(1,2)
    Poly(end,:)=[];
end

% find all the edges and number them
if isempty(PolyIn)
    vernum = size(Poly,1);
    Edge = [1:vernum;[2:vernum 1]]';
    PT=Poly;
else
    vernum1 = size(Poly,1);
    Edge = [1:vernum1;[2:vernum1 1]]';
    PT=Poly;
    for poi= 1:length(PolyIn)
        vernum2 = vernum1+size(PolyIn{poi},1);
        Edge = [Edge;[vernum1+1:vernum2;[vernum1+2:vernum2 vernum1+1]]'];
        PT = [PT;PolyIn{poi}];
        vernum1 = size(Edge,1);
    end    
end

% find a box containing all the polys and meshgrid the box
xymin=min(PT);
xymax=max(PT);
lx=xymax(1)-xymin(1);
ly=xymax(2)-xymin(2);
if isempty(dl)
    dl=max(lx,ly)/20;
end
nx=round(lx/dl);
ny=round(ly/dl);
dx=lx/nx;
dy=ly/ny;
x = xymin(1):dx:xymax(1);
y = xymin(2):dy:xymax(2);
[x0,y0]= meshgrid(x,y);
Nodes = [x0(:) y0(:)];

% find the nodes in the geometry
NodesInNum = find(inpolygon(Nodes(:,1),Nodes(:,2),Poly(:,1),Poly(:,2)));
if ~isempty(PolyIn)
    for poi = 1:length(PolyIn)
        NodesInNum0 = find(inpolygon(Nodes(:,1),Nodes(:,2),PolyIn{poi}(:,1),PolyIn{poi}(:,2)));
        NodesInNum = setdiff(NodesInNum,NodesInNum0);
    end
end

% delete the nodes very close to the edges, avoiding distortion element
% the control distance here is 0.2*dl
nodedel = zeros(length(NodesInNum),1);
for ii = 1:length(NodesInNum)
    cds = Nodes(NodesInNum(ii),:);
     % reference the function p_poly_dist to determine the distance from
     % nodes to edges
    d = p_poly_dist(cds(1),cds(2), Poly(:,1), Poly(:,2));
    if abs(d)<0.2*dl
        nodedel(ii)=1;
    end
end
for jj = 1:length(PolyIn);
    Polyjj = PolyIn{jj};
    for ii = 1:length(NodesInNum)
        cds = Nodes(NodesInNum(ii),:);
        d = p_poly_dist(cds(1),cds(2), Polyjj(:,1), Polyjj(:,2));
        if abs(d)<0.2*dl
            nodedel(ii)=1;
        end
    end
end
NodesInNum = setdiff(NodesInNum,NodesInNum(nodedel==1));

% divided all the edges to several segments
NodesIn = Nodes(NodesInNum,:);
NodesOnEdge=[];
for Edgenum = 1:length(Edge)
    vertex = PT(Edge(Edgenum,:),:);
    Le = norm(vertex(1,:)-vertex(2,:));
    num = ceil(Le/dl);
    xcds = vertex(1,1)+(vertex(2,1)-vertex(1,1))*(0:num)/num;
    ycds = vertex(1,2)+(vertex(2,2)-vertex(1,2))*(0:num)/num;
    NodesOnEdge = [NodesOnEdge;[xcds;ycds]'];
end

% combine the nodes and delaunay them
NodesAll = [NodesIn;NodesOnEdge];
NodesAll = unique(NodesAll,'rows');
tri = delaunay(NodesAll(:,1),NodesAll(:,2));

% determine if the tris are in the geometry,delete the outer one
triCenter = zeros(size(tri,1),2);
for triNum = 1:size(tri,1)
    triCenter(triNum,:) = sum(NodesAll(tri(triNum,:),:),1)/3;
end
triValid = find(inpolygon(triCenter(:,1),triCenter(:,2),Poly(:,1),Poly(:,2)));
if ~isempty(PolyIn)
    for poi = 1:length(PolyIn)
        triValid0 = find(inpolygon(triCenter(:,1),triCenter(:,2),PolyIn{poi}(:,1),PolyIn{poi}(:,2)));
        triValid = setdiff(triValid,triValid0);
    end
end 
TRI = tri(triValid,:);
PT = NodesAll;

% plot the TRI
figure
triplot(TRI,NodesAll(:,1),NodesAll(:,2));
%  triplot(TRI,NodesAll(:,1),NodesAll(:,2),'color','r','Marker','o','MarkerFaceColor','g');
axis equal 


function d = p_poly_dist(x, y, xv, yv) 
%*******************************************************************************
% function:	p_poly_dist
% Description:	distance from piont to polygon whose vertices are specified by the
%              vectors xv and yv
% Input:  
%    x - point's x coordinate
%    y - point's y coordinate
%    xv - vector of polygon vertices x coordinates
%    yv - vector of polygon vertices x coordinates
% Output: 
%    d - distance from point to polygon (defined as a minimal distance from 
%        point to any of polygon's ribs, positive if the point is outside the
%        polygon and negative otherwise)
% Routines: p_poly_dist.m
% Revision history:
%    7/9/2006  - case when all projections are outside of polygon ribs
%    23/5/2004 - created by Michael Yoshpe 
% Remarks:
%*******************************************************************************

    % If (xv,yv) is not closed, close it.
    xv = xv(:);
    yv = yv(:);
    Nv = length(xv);
    if ((xv(1) ~= xv(Nv)) || (yv(1) ~= yv(Nv)))
        xv = [xv ; xv(1)];
        yv = [yv ; yv(1)];
        Nv = Nv + 1;
    end

    % linear parameters of segments that connect the vertices
    A = -diff(yv);
    B =  diff(xv);
    C = yv(2:end).*xv(1:end-1) - xv(2:end).*yv(1:end-1);

    % find the projection of point (x,y) on each rib
    AB = 1./(A.^2 + B.^2);
    vv = (A*x+B*y+C);
    xp = x - (A.*AB).*vv;
    yp = y - (B.*AB).*vv;

    % find all cases where projected point is inside the segment
    idx_x = (((xp>=xv(1:end-1)) & (xp<=xv(2:end))) | ((xp>=xv(2:end)) & (xp<=xv(1:end-1))));
    idx_y = (((yp>=yv(1:end-1)) & (yp<=yv(2:end))) | ((yp>=yv(2:end)) & (yp<=yv(1:end-1))));
    idx = idx_x & idx_y;

    % distance from point (x,y) to the vertices
    dv = sqrt((xv(1:end-1)-x).^2 + (yv(1:end-1)-y).^2);

    if(~any(idx)) % all projections are outside of polygon ribs
       d = min(dv);
    else
       % distance from point (x,y) to the projection on ribs
       dp = sqrt((xp(idx)-x).^2 + (yp(idx)-y).^2);
       d = min(min(dv), min(dp));
    end

    if(inpolygon(x, y, xv, yv)) 
       d = -d;
    end
end

end