function pt=crvdirect(nodes, pt, dr)

% crvdirect: Check the direction of a curve and direct it.
% 
% Calling Sequences:
% 
%     pt=crvdirect(nodes, pt, dr)
% 
% INPUTS:
%
%      nodes - Coordinates of a set of 2D nodes.
%
%      pt  -  Indexes of the nodes.
%
%      dr - Direction of the curve. The default values is 1, 
%             which means the curve is anticlockwise. If dr=-1,
%             The curve will be direct to be clockwise.
%
% OUTPUT:
% 
%      pt - Indexes of the nodes after directing.
%  

[m, n]=size(nodes);
if n~=2 && m==2
    nodes=nodes';
end
if n~=2 && m~=2
    error('The nodes should be two-dimensional!');
end

if nargin==2
    dr=1;
end
cv=boundary(nodes(pt,1), nodes(pt,2), 0);
ec=[pt(cv(1:end-1)), pt(cv(2:end))];
ep=[pt(1:end-1), pt(2:end)];
de=min(DistanceMatrix(ec, ep), [], 2);
ii=find(de==0, 1);
if isempty(ii)
    pt=pt(end:-1:1);
end

% Make the direction to be clockwise
if dr==-1    
    pt=pt(end:-1:1);
end


