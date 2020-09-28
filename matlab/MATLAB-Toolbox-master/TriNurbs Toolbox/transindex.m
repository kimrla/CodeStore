function tri=transindex(tri,in)
% Find new indices in tri which are corresponding to all the nodes in order to corresponding to the points in the available domain.

% Input:
%   tri:Triangulations whose indices are corresponding to all the nodes.
%   in:Logical vector of one trimming surface in which 1 stands the point in the available domain and
%       0 out of the available domain.
% Output:
%   ntri:New indices of triangulations, which are corresponding to points
%       in the available domian instead of all the nodes. And triangulars whose
%       vertex doesn't belong to the domain are removed.



for j=1:3
    biaoji=[];
    n=length(tri);
    for i=1:n
        deter=in(tri(i,j));
        if (deter==1)
            tri(i,j)=sum(in(1:tri(i,j)));
        else
            biaoji=[biaoji,i];
        end
    end
    if (~isempty(biaoji))
        tri(biaoji,:)=[];
    end
end


end