function [lg,uv,st]=traceend_boundary(pts1,pts2,tol,tol34,varargin)
% Ending condition 1: When tracing to one of the 4 boundaries in the
% parameter domains of the 2 NURBS surfaces, if the par-coord is located 
% strictly at the boundary, then lg=1 and uv or st ={[],[]}; if the 
% par-coord is just near the boundary, then calculate the knot vector which
% is used in next ite3par.m.

% Input:
%   pts1,pts2: The 2 par-coords of the 2 NURBS surfaces, which is going to 
%       be checked. (Their corresponding phy-coord is the same).
%   tol: Tolerance, options.tolerance.
%   tol34: Tolerance to determine whether to transform from ite4par to
%       ite3par, options.tolerance34.
%   N: If transform to ite3par (meaning the tol34 has been satisfied), then
%       the refined knot vector. Default is 5. (linspace(a,b,N))

% Output:
%   lg: 0-Doesn't arrive to the boundary, maybe just near the boudnary;
%       1-Has been arrived to the boudnary of at least one of the surfaces..
%   uv,st: Refined knot vectors of the 2 surfaces respectively to use 
%       ite3par. If lg=1, then uv={[],[]} and st={[],[]}; If lg=0 and the 
%       point is near the boundary (satisfying the tol34), then refine the
%       knot interval of the 2 surfaces. And when use ite3par next, just
%       select one knot vector (u or v or s or t) as the iso-par-line.

if ~isempty(varargin)
    N=varargin{1};
else
    N=5;
end

lg=0;
pts=[pts1(:),pts2(:)];
uvst={[],[];[],[]};
for i=1:2
    if any(pts(:,i)<=tol) || any(abs(pts(:,i)-1)<=tol)
        % u or v is located at u=0 or 1 or v=0 or 1 strictly.
        lg=1;
        uv={[],[]};
        st={[],[]};
        return;
    else
        for j=1:2
            if pts(j,i)<=tol34
                uvst{j,i}=linspace(0,pts(j,i),N);
            elseif abs(pts(j,i)-1)<=tol34
                uvst{j,i}=linspace(pts(j,i),1,N);
            end
        end
    end
end

uv=uvst(1:2);
st=uvst(3:4);

end






