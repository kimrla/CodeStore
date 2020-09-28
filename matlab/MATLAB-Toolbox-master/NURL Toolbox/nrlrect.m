function curve = nrlrect(w,h)
% 
% NRLRECT: Construct NURL representation of a rectangular curve.
% 
% Calling Sequence:
% 
%   crv = nrlrect()
%   crv = nrlrect(size)
%   crv = nrlrect(width, height)
% 
% INPUT:
% 
%   size	: Size of the square (width = height).
% 
%   width	: Width of the rectangle (along x-axis).
% 
%   height	: Height of the rectangle (along y-axis).
%
% OUTPUT:
%
%   crv		: NURL curve, see nrlmake. 
%  
% 
% Description:
% 
%   Construct a rectangle or square in the x-y plane with the bottom
%   lhs corner at (0,0,0). If no rhs arguments provided the function
%   constructs a unit square.
%

if nargin < 1
   w = 1;
   h = 1;
end

if nargin < 2
   h = w;
end

coefs  = [0 w w 0 0;
               0 0  h  h 0;
               0 0  0  0 0;
               1 1  1  1 1];

knots  = [0 0.25 0.5 0.75 1];
intervals  = [0 0.25 0.5 0.75 1];

curve = nrlmake(coefs, knots, intervals, 1);

end

%% demo
%! crv = nrltform(nrbrect(2,1), vecrotz(deg2rad(35)));
%! nrbplot(crv,4);
%! axis equal
%! title('Construction and rotation of a rectangular curve.');
%! hold off





