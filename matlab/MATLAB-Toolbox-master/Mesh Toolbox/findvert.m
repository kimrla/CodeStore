function m=findvert(pnt,pnts,st1,st2)
% Find which line the point in parameter grids belongs to.
% input:
%   pnt: coordinates of the point in parameter domain
%   pnts: coordinates of each vertices of the parameter domain after
%         triangulation, see Triangulation.Points
%   st1: knots vector in x direction
%   st2: knots vector in y direction
% output:
%   m: ID of the vertices in one line, which is a 2-element vector

M=length(st1);
N=length(st2);
spx=st1(2)-st1(1);
spy=st2(2)-st2(1);

if (st1(1)==pnt(1))
       x1=0;x2=2;
elseif (st1(M)==pnt(1))
       x1=M-1;x2=M+1;
else 
    for i=1:M-1
        if (st1(i)<pnt(1) && st1(i+1)>pnt(1))
            x1=i;x2=i+1;
            break;
        elseif (st1(i)<pnt(1) && st1(i+1)==pnt(1))
            x1=i;x2=i+2;
            break;
        end
    end
end


if (st2(1)==pnt(2))
       y1=0;y2=2;
elseif (st2(N)==pnt(2))
       y1=N-1;y2=N+1;
else 
    for i=1:N-1
        if (st2(i)<pnt(2) && st2(i+1)>pnt(2))
            y1=i;y2=i+1;
            break;
        elseif (st2(i)<pnt(2) && st2(i+1)==pnt(2))
            y1=i;y2=i+2;
            break;
        end
    end
end

    %没有考虑当交点恰好位于纵横网格线交点的情况
if (x2-x1==2 && y2-y1==1)
    m(1)=find(pnts(:,1)==pnt(1) & pnts(:,2)<pnt(2) & pnts(:,2)>pnt(2)-spy);
    m(2)=find(pnts(:,1)==pnt(1) & pnts(:,2)<pnt(2)+spy & pnts(:,2)>pnt(2));
elseif (y2-y1==2 && x2-x1==1)
    m(1)=find(pnts(:,2)==pnt(2) & pnts(:,1)<pnt(1) & pnts(:,1)>pnt(1)-spx);
    m(2)=find(pnts(:,2)==pnt(2) & pnts(:,1)<pnt(1)+spx & pnts(:,1)>pnt(1));
end
end




