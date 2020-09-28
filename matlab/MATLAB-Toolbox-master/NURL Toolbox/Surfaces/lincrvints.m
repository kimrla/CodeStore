function [u, d, pp]=lincrvints(crv, lin, u, it)

% Solve a intersection point of two curves using Newton-Raphson's method
%
%  Inputs: 
% 
%    crv, lin - the curve and the line to be used
% 
%    u - an initial parametric point of the two curves
%
%  Output: 
% 
%    u - the intersection parametric point
% 
%    d - the minimum distance of the two curves
% 
%    pp - the coordinates of the intersection point of the two curves
%

if nargin==3
    it=5;
end

% Solve the intersectiont point using Newton-Raphson's method
F=zeros(2,1); dF=zeros(2);
for i=1:it
    [p1, dp1, hessm1] = nrldeval(crv, u(1)); 
    [p2, dp2] = nrldeval(lin, u(2)); 
    F(1)=2*sum((p1-p2).*dp1); 
    F(2)=-2*sum((p1-p2).*dp2); 
    dF(1,1)=2*sum(dp1.^2)+sum((p1-p2).*(hessm1));
    dF(1,2)=-2*sum(dp1.*dp2);
    dF(2,1)=dF(1,2);
    dF(2,2)=2*sum(dp2.^2);
    u=u-dF\F;
    for j=1:2 
        if u(j)<0 
            u(j)=0; 
        elseif u(j)>1 
            u(j)=1; 
        end
    end
end
d=sqrt(sum((p1-p2).^2)); 
pp=[p1,p2];
        







