function [q, p]=apntsints(pnts1, pnts2, tol)

% Nearest distances of two sets of points
%
% Calling Sequences:
% 
%     [p, q]=asrfcrvints(pnts1, pnts2, tol)
% 
% INPUTS:
%
%     pts1 - points of intersections on the surface
%     pts2 - points of intersections on the curve
%     tol - tolerance of length
% 
% OUTPUT:
% 
%    q - index of intersection points on pnts1
%    p - index of intersection points on pnts2
% 
% See also:
%    srfsintersct
%  

k=length(pnts2(1,:));
dm=DistanceMatrix(pnts2(:,:)', pnts1(:,:)'); 
pp=dm<tol;
qq=zeros(k, 1);
for i=1 : k 
    qi=find(pp(i, :)); 
    if ~isempty(qi) 
        di=DistanceMatrix(pnts2(:,i)', pnts1(:,qi)'); 
        [~, id]=min(di); 
        qq(i)=qi(id); 
    end
end
id=qq~=0; 
p=find(id); q=qq(id); 







