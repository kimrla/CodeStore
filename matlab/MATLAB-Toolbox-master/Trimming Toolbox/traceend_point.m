function [lg,uv]=traceend_point(pts1,pnts1,ptslist,pntslist,tol,tol34,varargin)
% Ending condition 2: Check if the 2 points pnts1 and the ending point in 
% intersection point list ptslist are coincident, or if pnts1 has 'entered'
% into the pntslist (meaning the intersection part has been connected as a
% closed circle or one part). If the pnts1 is Not coincident, but very near
% the ending point, then refine the knot vector to use ite3par next. Notice
% that all the par-coords should be selected in the same surface (srf1 or
% srf2), and don't change the par-domain.

% Input:
%   pts1,pnts1: Par/phy-coords of the checking piont, where the par-coord
%       is in srf1 or srf2.
%   ptslist,pntslist: Par/phy intersection point list, where the par-coords
%       are in srf1 or srf2. And the FIRST point in the list is the
%       checking ending point. size(list)=[numpoint,2 or 3]
%   tol: Tolerance, options.tolerance.
%   tol34: Tolerance to determine whether to transform from ite4par to
%       ite3par, options.tolerance34.
%   N: If transform to ite3par (meaning the tol34 has been satisfied), then
%       the refined knot vector. Default is 5. (linspace(a,b,N))

% Output:
%   lg: 0-Doesn't arrive to the ending point, maybe just near the ending;
%       1-Has been coincident to the ending point, or 'enter' the ptslist.
%   uv: Refined knot vectors of ONE of the 2 surfaces to use 
%       ite3par. If lg=1, then uv={[],[]} and st={[],[]}; If lg=0 and the 
%       point is near the ending point (satisfying the tol34), then refine
%       knot interval of ONE of the 2 surfaces (srf1 or srf2, depending on 
%       the input pts1 and ptslist). And when use ite3par next, just
%       select one knot vector (u or v) as the iso-par-line. If the point
%       pts1 'enter' the list, then uv=d, which is the distance.

if ~isempty(varargin)
    N=varargin{1};
else
    N=5;
end

pts2=ptslist(1,:);
pnts2=pntslist(1,:); % The ending point of the intersection list
pnts3=pntslist(2,:); % The 2nd point in the intersection list

lg=0;
uv={[],[]};

if norm(pnts1-pnts2)<=tol
    lg=1;           
    return;
else
    tem1=dot(pnts3-pnts2,pnts1-pnts2);
    if tem1>0 % The last point has 'enter' the first and second points.
        % Fit the intersection curve, and project the point into the curve.
        crv=bspinterpcrv(pntslist',2);% Default order of the curve is 2.
        [~,~,d2]=projpt2crv(crv,pts1,pnts1);
        if d2<=tol*10% NOTICE: The colerance can be adjusted.
            lg=1;
            uv=d2;
            return;
        end
    else
        if norm(pnts1-pnts2)<=tol34% tol34 is related to ite4par to ite3par.
            uv{1}=sort(linspace(pts1(1),pts2(1),N));
            uv{2}=sort(linspace(pts1(2),pts2(2),N));
        end
    end
end

end



