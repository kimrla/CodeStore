function n=unitnvector(srf,pts)
% Calculate the unit normal vector of a point on the NURBS surface. If the
% tangent vector cannot be obtained by the derivative, then use difference
% quotient to replace it automatically. 

% Input:
%   srf: NURBS/Bezier surface.
%   pts: Par-coords of the given point on the surface.

% Output:
%   n: Unit normal vector of the point pts.

f=@(knt) (nrbeval(srf,knt));
dsrf=nrbderiv(srf);
lg=isderiv(dsrf);

if lg==1
    dF=discretederiv(f,pts); % size(dF)=[3,2]
else
    [~,dF]=nrbdeval(srf,dsrf,pts);   
    dF=[dF{:}];
end
n=vecnorm(cross(dF(:,1),dF(:,2)));

end