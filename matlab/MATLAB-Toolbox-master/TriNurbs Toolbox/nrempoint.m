function [points,logvec,id]=nrempoint(pnts,tol)

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
bas=1;
id=[];

% Firstly,check distances between every 2 adjacent points along the curve,if the
% distance<tol,then remove the latter point(record the index of the latter
% point)
while (bas<n)
    k=bas;
    for i=k+1:n
        d=norm(pnts(i,:)-pnts(k,:));
        if (d<tol)
           id(end+1)=i;
           bas=i;
        else
            bas=bas+1;
            break;
        end
    end
end

temp=ones(1,n);
if (~isempty(id))
    temp(id)=0;
    logvec=logical(temp);
    points=pnts(logvec,:);
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
% % Plot results
% figure; hold on; 
% tnrbplot(tsrf);
% tnrbplot(tcrv);
% plot(pnts1(:,1), pnts1(:,2), 'ro');
% plot(pnts2(:,1), pnts2(:,2), 'r*');
% axis equal;
% figure; hold on;
% triplot(tsrf.delaunay, tsrf.nodes(:,1), tsrf.nodes(:,2)); 
% plot(pts1(:,1), pts1(:,2), 'k.', 'MarkerSize', 13); 
% title('Parametric mesh of the surface'); 
% axis equal;
% 
% % Firstly remove the near points just along the curve,then remove the
% % near points between nodes of curve and the surface.
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
% tempoints(remindex,:)=[];%remove redundant points in physical domain of surface
% tempoints=[tempoints;newpnts];%newpnts:remove redundant points in physical domain of curve
% temnodes=tsrf.nodes;
% temnodes(remindex,:)=[];%%remove redundant points in parameter domain of surface
% temnodes=[temnodes;pts1(logvec,:)];%pts1(logvec,:):remove redundant points in parameter domain of curve

% tridel=delaunayTriangulation(temnodes);
% temp2=tridel.ConnectivityList;
% figure;
% trisurf(temp2,tempoints(:,1),tempoints(:,2),tempoints(:,3));
% axis equal;view(2);
% figure;
% triplot(temp2,temnodes(:,1),temnodes(:,2));
% axis equal;view(2);


            
            
            
            
            
            
            
            