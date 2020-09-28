function [npnts,nbndin]=delaunaydeal(nrb,pts,bndin,pnts,oldbndin)
% This function is used to deal with the new points when the MATLAB
% function delaunayTriangulation inserts or includes new points in the
% boundary.Only deal with a continuous curve,including the
% self-intersection condition.

% Input:
%   nrb: Nurbs surface structure.
%   pts: Parameter coordinates after the execution of 
%       tr=delaunayTriangulaiton,see tr.Points.
%   bndin: Indices of points in the boundary corresponding to npts. Each
%       row represents a boudnary edge, but the order may not be sorted, see tr.Constraints.
%   pnts: Physical coordinates of the points before the execution of
%       tr=delaunayTriangulation, which needs to be dealt with.
%   oldbndin: Indices of points in the boundary corresponding to the
%       original parametric points (or pnts) before the function
%       delaunayTriangulation inserts or includes new points.

% Output:
%   npnts: New physical coordinates after dealing with the inserting
%       points.
%   nbndin: New Indices of boundary points, which is rearranged into a sorted
%       sequence to form a boudnary.

m=length(pts);
n=length(pnts);
tembndin=bndin;
nbndin=[];
numedge=length(oldbndin);

for i=1:numedge
    % Deal with each edge in oldbndin
    startbnd=oldbndin(i,1);
    endbnd=oldbndin(i,2);
    while (true)
        % Find point connecting with startbnd
        [u,v]=find(tembndin==startbnd);
        k=2./v;
        anotherstartbnd=[];
        for ij=1:length(u)
            anotherstartbnd=[anotherstartbnd,tembndin(u(ij),k(ij))];
        end
        % Cancel the condition of connection between the points with startbnd
        paid=anotherstartbnd==startbnd;
        u(paid)=[];v(paid)=[];anotherstartbnd(paid)=[];
        k=2./v;
        numu=length(u);
        % Condition of no self-intersection
        if (numu==1)
            % The last new inducing point in the edge
            if (anotherstartbnd==endbnd)
                nbndin=[nbndin;startbnd,endbnd];
                tembndin(u,:)=[];
                break;
            % Insert the new inducing point in the edge
            else
                nbndin=[nbndin;startbnd,anotherstartbnd];
                startbnd=anotherstartbnd;
                tembndin(u,:)=[];
            end
        % Condition of the self-intersection point.Determine whether the point is located in the edge    
        elseif (numu>1)
            linepoint1=pts(startbnd,:);
            linepoint2=pts(endbnd,:);
            for j=1:numu
                trypoint=pts(anotherstartbnd(j),:);
                z=(trypoint(2)-linepoint1(2))*(linepoint2(1)-linepoint1(1))-(linepoint2(2)-linepoint1(2))*(trypoint(1)-linepoint1(1));
                z=abs(z);
                if (z<1e-6)
                    nbndin=[nbndin;startbnd,anotherstartbnd(j)];
                    startbnd=anotherstartbnd(j);
                    tembndin(u(j),:)=[];
                    break;
                end
            end
            if (startbnd==endbnd)
                break;
            end
        % The endding point of the curve and the curve is not a closed
        % loop,which means the last point doesn't connect with the first one
        else
            nbndin=[nbndin;startbnd,endbnd];
            break;
        end
    end
end
                    
% Deal with the inserting points
% The boundary just consists of new points which exists in the original surface
if (m==n) 
    npnts=pnts;
%the boundary inserts new points into the boundary
else 
    yaoqiude=pts(n+1:end,:)';
    p=nrbeval(nrb,yaoqiude);
    p=p(:,:)';
    npnts=[pnts;p];    
end
                
end

% %demo
% n=10;
% x=1:n;
% y=1:n;
% [X,Y]=meshgrid(x,y);
% 
% bnd=[12,33,32,53,54,57,15,64];
% bnd=[bnd(1:end);bnd(2:end),12]';
% 
% tr=delaunayTriangulation(X(:),Y(:),bnd);
% 
% pnt=tr.Points;
% tri=tr.ConnectivityList;
% cnt=tr.Constraints;
% 
% 
% m=length(cnt);
% figure;
% hold on;
% triplot(tri,pnt(:,1),pnt(:,2));
% plot(pnt(bnd,1),pnt(bnd,2),'ro','Markersize',8);
% for i=1:m
%     plot(pnt(cnt(i,:),1),pnt(cnt(i,:),2),'dg','LInewidth',2,'Markersize',5);
% end
% 
% nbndin=delaunaydeal(nrb,pnt,cnt,bnd,bnd);%need to add the nurbs model

