function [ed1, ed2, pts1, pts2, pnts1, pnts2]=tnrbinterline(tsrf, tcrv, p2t1, p1, p2)

% tnrbinterline: Get the intersection points of a tri-nurbs surface with a curve
% 
% Calling Sequences:
% 
%       [ed1, ed2, pts1, pts2, pnts1, pnts2]=tnrbinterline(tsrf, tcrv, p2t1, p1, p2)
% 
% INPUTS:
% 
%       tsrf, tcrv - Triangular representation of nurbs surface and curve.
% 
%       p2t1 - The relations from points to triangles of tri-nurbs 
%                  surface 1. See also tnrbpts2tri.
% 
%       p1, p2 - Indexes of nearest points between the surface and the curve.
%                    See also nearpnts.
%
% OUTPUT: 
%
%      ed1  -  Edges of the triangles of tri-nurbs surface 1 that
%                   intersected with the curve.
%
%      ed2   -  Line segments of tri-nurbs curve that intersecte with 
%                  tri-nurbs surface 1.
%
%       pts1, pts2 - Parametric intersection points of the surface and the curve.
%
%       pnts1, pnts2 - Intersection points of the surface and the curve.
%

tol=max([tsrf.seeds(1), tcrv.seeds(1)]);
nd=length(p1); t=0;
for k=1:nd 
    edges=tnrbpt2edges(tsrf, p2t1, p1(k));
    m=length(edges);
    for i=1:m
        ei1=[p1(k), edges(i)];
        if p2(k)==1
            ei2=[p2(k), p2(k)+1; p2(k)+1, p2(k)+2];
        elseif p2(k)==2
            ei2=[p2(k)-1, p2(k); p2(k), p2(k)+1; p2(k)+1, p2(k)+2];
        elseif p2(k)==tcrv.numbers(1)
            ei2=[p2(k)-2, p2(k)-1; p2(k)-1, p2(k)];
       elseif p2(k)==tcrv.numbers(1)-1
           ei2=[p2(k)-2, p2(k)-1; p2(k)-1, p2(k); p2(k), p2(k)+1];
        else
            ei2=[p2(k)-2, p2(k)-1; p2(k)-1, p2(k); p2(k), p2(k)+1; p2(k)+1, p2(k)+2];
        end
        L1=tsrf.points(ei1, :); 
        for j=1:size(ei2,1)
            L2=tcrv.points(ei2(j,:), :); 
            [tpnts1, tpnts2, x, d]=interline2line(L1, L2);
            p= x(1,:)>-eps & x(1,:)<1+eps; 
            p=p & d<tol*1e-6;
            p=p & x(2,:)>-eps & x(2,:)<1+eps; 
            x=x(:,p); 
            tpnts1=tpnts1(p,:); 
            tpnts2=tpnts2(p,:); 
            line1=tsrf.nodes(ei1, :); 
            line2=tcrv.nodes(ei2(j,:), :); 
            if ~isempty(x)
                tpts1x=(1-x(1,:))*line1(1,1)+x(1,:)*line1(2,1); 
                tpts1y=(1-x(1,:))*line1(1,2)+x(1,:)*line1(2,2); 
                tpts1=[tpts1x(:), tpts1y(:)]; 
                tpts2=(1-x(2,:))*line2(1)+x(2,:)*line2(2); 
                ed1(t+1,:)=ei1;
                ed2(t+1,:)=ei2(j,:);
                pts1(t+1,:)=tpts1(1,:);
                pts2(t+1,:)=tpts2(1,1);
                pnts1(t+1,:)=tpnts1(1,:);
                pnts2(t+1,:)=tpnts2(1,:);
                t=t+1;
            end
        end
    end
end
ed1=sort(ed1, 2);
[ed1, id]=RemDuplicate(ed1);
ed2=ed2(id,:);
pts1=pts1(id,:);
pts2=pts2(id,:);
pnts1=pnts1(id,:);
pnts2=pnts2(id,:);
[pts2, id]=sort(pts2);
ed1=ed1(id,:);
ed2=ed2(id,:);
pts1=pts1(id,:);
pnts1=pnts1(id,:);
pnts2=pnts2(id,:);

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
% 
% figure; hold on;
% triplot(tsrf.delaunay, tsrf.nodes(:,1), tsrf.nodes(:,2)); 
% plot(pts1(:,1), pts1(:,2), 'k.', 'MarkerSize', 13); 
% title('Parametric mesh of the surface'); 
% axis equal;





