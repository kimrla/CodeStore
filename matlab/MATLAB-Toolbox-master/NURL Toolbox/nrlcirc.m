function curve = nrlcirc(radius,center,sang,eang)
% 
% NRLCIRC: Construct a circular arc.
% 
% Calling Sequence:
% 
%   crv = nrlcirc()
%   crv = nrlcirc(radius)
%   crv = nrlcirc(radius,center)
%   crv = nrlcirc(radius,center,sang,eang)
% 
% INPUT:
% 
%   radius	: Radius of the circle, default 1.0
% 
%   center	: Center of the circle, default (0,0,0)
% 
%   sang	: Start angle, default 0 radians (0 degrees)
% 
%   eang	: End angle, default 2*pi radians (360 degrees)
% 
% OUTPUT:
%
%   crv		: NURL curve for a circular arc.
% 
% Description:
% 
%   Constructs NURL data structure for a circular arc in the x-y plane. If
%   no rhs arguments are supplied a unit circle with center (0.0,0.0) is
%   constructed. 
% 

if nargin < 1
  radius = 1;
end

if nargin < 2
  center = [];
end
   
if nargin < 4
  sang = 0;
  eang = 2*pi;
end

curve = nrlellip(radius, radius, [0, 0], sang, eang);

% vectrans arc if necessary
if ~isempty(center) 
  xx = vectrans(center);
  curve.coefs(1:3, :)=curve.coefs(1:3, :).*repmat(curve.coefs(4, :), 3, 1);
  curve.coefs = xx*curve.coefs;
  curve.coefs(1:3, :)=curve.coefs(1:3, :)./repmat(curve.coefs(4, :), 3, 1);
end

end

%% demo
%  for r = 1:9
%     crv = nrlcirc(r,[],deg2rad(45),deg2rad(315));
%     nrlplot(crv,50);
%     hold on;
%  end
%  hold off;
%  axis equal;
%  title('NURL construction of several 2D arcs.');
