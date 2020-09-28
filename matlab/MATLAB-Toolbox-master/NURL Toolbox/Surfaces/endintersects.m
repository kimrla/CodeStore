function [x, pnts1, pnts2, dt]=endintersects(srf1, srf2, x, pnts1, pnts2, dt, p)

% Add boundary intersection points of two surfaces to intersects
% 
% Calling Sequences:
% 
%     [x, pnts1, pnts2, dt]=endintersects(srf1, srf2, x, pnts1, pnts2, dt, p)
% 
% INPUTS:
%
%     srf1, srf2 - nurls surfaces 
%     x - parametric points of intersections
%     pnts1 - points of intersections on surface 1
%     pnts2 - points of intersections on surface 2
%     dt - distance of the two nearest points
%     p - index of boundary points
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     pnts1 - points of intersections on surface 1
%     pnts2 - points of intersections on surface 2
%     dt - distance of the two nearest points
% 
%  See also:
%       srfsinterscts
% 

n=length(p);
for k=1:n
    ii=mod(p(k), 4);
    jj=fix(p(k)/4);
    if ii==0
        ii=4;
    else
        jj=jj+1;
    end
    if ii==1
        crv=nrlsrf2crv(srf1, 0, []); 
        xi=x(2:4,jj);
        [xi, pntsi1, pntsi2, dti]=srfcrvintersect(srf2, crv, xi, 3);
        xi=[0; xi];
    elseif ii==2
        crv=nrlsrf2crv(srf1, [], 0); 
        xi=x([1, 3, 4], jj);
        [xi, pntsi1, pntsi2, dti]=srfcrvintersect(srf2, crv, xi, 3);
        xi=[xi(1); 0; xi(2:3)];
    elseif ii==3
        crv=nrlsrf2crv(srf2, 0, []); 
        xi=x([1,2,4],jj);
        [xi, pntsi1, pntsi2, dti]=srfcrvintersect(srf1, crv, xi, 3);
        xi=[xi(1:2); 0; xi(3)]; 
    else
        crv=nrlsrf2crv(srf2, [], 0); 
        xi=x(1:3,jj);
        [xi, pntsi1, pntsi2, dti]=srfcrvintersect(srf1, crv, xi, 3);
        xi=[xi(1:3); 0]; 
    end
    x(:,jj)=xi(:);
    pnts1(:,jj)=pntsi1(:);
    pnts2(:,jj)=pntsi2(:);
    dt(jj)=dti;
end

