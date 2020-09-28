function m=iftri(pnt1,pnt2,st1,st2,tri)
% Determine whether the adjacent 2 points in the parameter grids belongs to
%   the same one triangular element
% Input:
%   tri: triangulation structure, see Triangulation or delaunayTriangulation
%   pnt1,pnt2: coordinates of the 2 adjacent points in parameter domain
%   st1: knots vector in x direction of the surface parameter domain
%   st2: knots vector in y direction of the surface parameter domain
% output:
%   m: Results of whether the 2 points belong to the same triangular, if
%       yes, then return the ID of the triangular(See Triangulation.ConnectivityList),
%       if not, return a 1*4 vector which contains 4 IDs of the grid
%       containing the 2 adjacent points

spx=st1(2)-st1(1);
spy=st2(2)-st2(1);
pnts=tri.Points;

% �����ǽ���ǡ�ô���������غϵ����
%�Ա�������Խ��߱��н���
if (abs(abs(pnt1(1)-pnt2(1))-spx)<1e-6 || abs(abs(pnt1(2)-pnt2(2))-spy)<1e-6)
    m1=findvert(pnt1,pnts,st1,st2);
    m2=findvert(pnt2,pnts,st1,st2);
    m=sort([m1,m2]);
else
    %�ڱ����
    [m1,linid1]=findtri(tri,pnt1,st1,st2);
    [m2,linid2]=findtri(tri,pnt2,st1,st2);
%     n1=length(m1);
%     n2=length(m2);
    m=sort([m1,m2]);
    m_=[unique(m),-1];
    %�ڱ��ҶԽ������н���
    if (length(m)~=length(m_))
        m=sort([linid1,linid2]);
    else
        %�ڱ��ҶԽ������޽��㣬������������ͬһ�������Σ�ֻ���ظ������ε�ID
        for i=1:length(m_)
            if (m(i)~=m_(i))
                m=m(i);
                break;
            end
        end
    end
end
end
        
           
                
        










