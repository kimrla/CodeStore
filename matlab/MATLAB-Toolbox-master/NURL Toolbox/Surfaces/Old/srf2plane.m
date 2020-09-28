function [pnts, qq]=srf2plane(pnts1, nor1, nor2, D10, D20, tol, ii)

% Get intersection of two surfaces by tangent planes
% 
% Calling Sequences:
%
%     [pnts, qq]=srf2plane(pnts1, nor1, nor2, D10, D20, tol, ii)
%
% INPUTS:
%
%      pnts1  - points of surfaces 1
%      nor1  -  normal vectors of surfaces 1
%      nor2  -  normal vectors of surfaces 2
%      D10, D20 - constants of plane equations of the two surfaces
%      tol - tolerace
%      ii - the coordinates that will not change (1, 2, or 3)
%
% OUTPUT:
% 
%      pnts  - points of or intersection
%      qq - index of bad intersection points
%

jj=1:3; jj(ii)=[];
D1=D10-nor1(ii,:).*pnts1(ii,:);
D2=D20-nor2(ii,:).*pnts1(ii,:);
pnts=pnts1;
[~, n]=size(pnts);
for i=1:n
    M=[nor1(jj(1),i), nor1(jj(2),i); nor2(jj(1),i), nor2(jj(2),i)]; 
    D=[D1(i); D2(i)];
    pnts([jj(1), jj(2)],i)=M\D;
end

qq=zeros(n, 1);
for i=1:n
    if norm(pnts(:,i)-pnts1(:,i))>tol
        qq(i)=i;
    end
end
qq=find(qq~=0);





