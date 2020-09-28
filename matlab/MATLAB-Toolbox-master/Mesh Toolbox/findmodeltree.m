function id=findmodeltree(modeltree,n,sign)
% Use the tree.node to find the corresponding node. sign can be
% 'depth','parent','current','child' and n is the corresponding index.

% id: Index of the corresponding modeltree: modeltree{id}.
% modeltree: All the nodes of the tree.
num=length(modeltree);     
id=[];
switch sign
    case 'depth'
        detid=1;% n: Depth of the searched node
    case 'parent'
        detid=2;% n: Parent node's index
    case 'current'
        detid=3;% n: Index of the current node
    case 'child'
        detid=4;% n: Number of children nodes
    otherwise
        error('The sign is wrong!');
end

for i=1:num
    node=modeltree{i}.node;
    if n==node(detid)
        id=[id,i];
    end
end

end