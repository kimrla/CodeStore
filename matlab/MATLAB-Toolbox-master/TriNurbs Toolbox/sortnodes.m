function num=sortnodes(pts, num)

% sortnodes: Sort the sequences of a few 2D nodes
% 
% Calling Sequences:
% 
%       num=sortnodes(pts, num)
% 
% INPUTS:
% 
%       pts - Parametric points of a surface or plane nodes.       
% 
%       num - Indexes of the node in pts.
%
% OUTPUT: 
% 
%       num - Indexes of the node in pts after sorting.
%
% Discription:
%
%      This routine only applies to a few nodes that have 
%      one coordinates monotonous. For example, nodes 
%      within a triangle.
%

[~, sv]=sort(pts(num,2));
num=num(sv);
for i=1:length(num)-2
    du1=vecnorm(pts(num(i+1),:)-pts(num(i),:));
    du2=vecnorm(pts(num(i+2),:)-pts(num(i+1),:));
    dr=dot(du1, du2);
    if dr<0
        [~, su]=sort(pts(num,1));
        num=num(su);
    end
end





