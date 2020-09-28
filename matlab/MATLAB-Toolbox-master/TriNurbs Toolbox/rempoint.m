function [points,logvec,id]=rempoint(pnts,tol)

% Remove the points whose distances are very small in a series of points.
% The input must be a series of points' coordinates that the points have
% been sorted and connected, see sortnodes, tnrbintersort and
% tnrbinterconnct.

% Input:
%   pnts: Sorted points, which is N*3 to represent the points' coordinates
%   tol: Checking standard
% Output:
%   points: Final points' coordinates after removing redundant points
%   logvec: 1*N logical vector, in which 1 means the corresponding point is
%       remained and 0 means the point is removed.
%   id: Index of the points which are removed. The order corresponds to the
%       input points.(1*M vector)

n=length(pnts);
k=1;
adjust=0;
id=-1;

for i=1:2:n-2
%     i=i-adjust;
    % distance of the 3 ponits with each other
    d1=norm(pnts(i-adjust,:)-pnts(i+1,:));
    d2=norm(pnts(i+2,:)-pnts(i+1,:));
    d3=norm(pnts(i-adjust,:)-pnts(i+2,:));
    [m,man]=sort([d1,d2,d3]);
    % need to remove 1 or 2 points
    if (m(1)<tol)
        % distance from 1 to 3 is the largest
        if (man(3)==3)
            id(k)=i+1;
            adjust=0;
            k=k+1;
            % continue to remove the second point
            if (m(2)<tol)
                id(k)=i+2;
                k=k+1;
                adjust=2;
            end
        else
            id(k)=i+2;
            adjust=1;
            k=k+1;
            if (m(2)<tol)
                id(k)=i+1;
                k=k+1;
                adjust=2;
            end
        end
    end
end

temp=ones(1,n);
% if (id==-1)
%     error('There is on point to be removed');
% else
if (k~=1)
    temp(id)=0;
    logvec=logical(temp);
    points=pnts(logvec,:);
end
% whether there are odd or even points, if there are even points, the last
% point need to be discussed
if (rem(n,2)==0)
    temp1=norm(points(end,:)-pnts(end,:));
    if (temp1<tol)
        id(k)=n;
        k=k+1;
        logvec(end)=0;
    else
        points(end+1,:)=pnts(end,:);
    end
end
    
%% demo
% % The mesh seed length (h0)
% h0=0.5;
% 
% % Create a plane square and a plane cuve
% lin1=nrbline([0,1], [9,1]);
% lin2=nrbline([0,6], [9,6]);
% srf=nrbruled(lin1, lin2);
% crv=nrbtestcrv;
% 
% % Transform a nurbs surface into triangular representation
% tsrf=nrb2tri(srf, h0);
% tcrv=nrb2tri(crv, h0);
% 
% % The nearest points from the surface to the curve
% tol=max([tsrf.seeds(1), tcrv.seeds(1)]);
% [p1, p2, d]=nearpnts(tsrf.points, tcrv.points, tol);
% 
% % Get the relations from points to triangles of tri-nurbs
% p2t1=tnrbpts2tri(tsrf);
% 
% % Get the intersection points of a tri-nurbs surface with a curve
% [ed1, ed2, pts1, pts2, pnts1, pnts2]=tnrbinterline(tsrf, tcrv, p2t1, p1, p2);       
% 
% % Plot results物理域
% figure; hold on; 
% tnrbplot(tsrf);
% tnrbplot(tcrv);
% plot(pnts1(:,1), pnts1(:,2), 'ro');
% plot(pnts2(:,1), pnts2(:,2), 'r*');
% axis equal;
% %参数域
% figure; hold on;
% triplot(tsrf.delaunay, tsrf.nodes(:,1), tsrf.nodes(:,2)); 
% plot(pts1(:,1), pts1(:,2), 'k.', 'MarkerSize', 13); 
% title('Parametric mesh of the surface'); 
% axis equal;
% 
% %先不考虑曲线的分段点,即先不反求
% tolerance=1/4*tol;
% [newpnts,logvec,id]=rempoint(pnts1,tolerance);
% [newp1,newp2,newd]=nearpnts(newpnts,tsrf.points);
% nn=length(newp1);
% remindex=-2;
% kk=1;
% for i=1:nn
%     if (newd(i)<tolerance)
%         remindex(kk)=newp2(i);
%         kk=kk+1;
%     end
% end

% tempoints=tsrf.points;
% tempoints(remindex,:)=[];%曲面物理域删除多余结点后
% tempoints=[tempoints;newpnts];%newpnts：交线物理域删除多余节点后
% temnodes=tsrf.nodes;
% temnodes(remindex,:)=[];%曲面参数域删除多余节点后
% temnodes=[temnodes;pts1(logvec,:)];%pts1(logvec,:)：交线参数域删除多余节点后

% tridel=delaunayTriangulation(temnodes);
% temp2=tridel.ConnectivityList;
% figure;
% trisurf(temp2,tempoints(:,1),tempoints(:,2),tempoints(:,3));
% axis equal;view(2);
% figure;
% triplot(temp2,temnodes(:,1),temnodes(:,2));
% axis equal;view(2);
% 
% 
% 
% 

            
            
            
            
            
            
            
            