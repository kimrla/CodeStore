function curve = nrlline(p1,p2)
% 
% NRLLINE: Construct a straight line.
% 
% Calling Sequence:
% 
%   crv = nrlline()
%   crv = nrlline(p1,p2)
% 
% INPUT:
% 
% p1		: 2D or 3D cartesian coordinate of the start point.
% 
% p2            : 2D or 3D cartesian coordinate of the end point.
%
% OUTPUT:
% 
% crv		: NURL curve for a straight line.
% 
% Description:
% 
%   Constructs NURL data structure for a straight line. If no rhs 
%   coordinates are included the function returns a unit straight
%   line along the x-axis.
%

coefs = [zeros(3,2); ones(1,2)];

if nargin < 2
  coefs(1,2) = 1.0;  
else
  coefs(1:length(p1),1) = p1(:);    
  coefs(1:length(p2),2) = p2(:);
end

curve = nrlmake(coefs, [0 1], [0, 1], 1);

end

%% demo
%  crv = nrlline([0.0 0.0 0.0]',[5.0 4.0 2.0]');
%  nrlplot(crv,2);
%  grid on;
%  title('3D straight line.');
%  hold off
