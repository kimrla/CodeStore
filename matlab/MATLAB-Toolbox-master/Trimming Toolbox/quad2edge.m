function qd2edge=quad2edge(qnrb,k)
% Find relationship between quad in qnrb.quad and edges in qnrb.qedges.

% Input:
%   qnrb:Quad-nurbs structure
%   k:The kth quad which needs to be searched

% Output:
%   qd2edge:Four edge indices corresponding to qnrb.qedges, whcih is a 1*4
%       array.

quad=qnrb.quad(k,:);
edges=qnrb.qedges;
edges=sort(edges,2);
% The four vertices of a quad element should be sorted in specific order,
% eg. the 2 vertices cannot be connected into a diagonal.

for j=1:3
    edge(j,:)=sort(quad([j,j+1]));
end
edge(j+1,:)=sort(quad([1,4]));

for j=1:4
    qd2edge(j)=find(edges(:,1)==edge(j,1) & edges(:,2)==edge(j,2));
end

end


    

