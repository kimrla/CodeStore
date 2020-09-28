function [uu, pp, d]=crvsints(crv1, crv2, dist, uu, du)

% Solve points of intersection of two curves
% 
%  Inputs: 
%
%    crv1, crv2 - the two curves to be used
% 
%    dist - tolerace
%
%     uu - the approximate intersection parametric points
%
%     du - the approximate minimum distances of the two curves
%
%  Output: 
%
%     uu - the intersection parametric points
%     pp - the coordinates of the intersection points of the two curves
%     d - the minimum distances of the two curves
% 

% Solve points of intersection using Newton-Raphson's method
m=length(du); d=zeros(m,1); pp=zeros(3,2,m); 
for i=1:m
    [uu(i,:), d(i), pp(:,:,i)]=intsctpnt(crv1, crv2, uu(i,:)');
end
id=d/dist<1e-6;
d=d(id); uu=uu(id,:); pp=pp(:,:,id);

[uu, dc]=RemDuplicate(uu);
pp=pp(:,:,dc); d=d(dc);






