function surf = nrlcylind(height,radius,center,sang,eang)
% 
% NRLCYLIND: Construct a cylinder or cylindrical patch.
% 
% Calling Sequence:
% 
%   srf = nrlcylind()
%   srf = nrlcylind(height)
%   srf = nrlcylind(height,radius)
%   srf = nrlcylind(height,radius,center)
%   srf = nrlcylind(height,radius,center,sang,eang)
% 
% INPUT:
% 
%   height	: Height of the cylinder along the axis, default 1.0
% 
%   radius	: Radius of the cylinder, default 1.0
% 
%   center	: Center of the cylinder, default (0,0,0)
% 
%   sang	: Start angle relative to the origin, default 0.
% 
%   eang	: End angle relative to the origin, default 2*pi.
%
% OUTPUT: 
%
%   srf     : cylindrical surface patch 
% 
% Description:
% 
%   Construct a cylinder or cylindrical patch by extruding a circular arc.
%


if nargin < 1
  height = 1;
end

if nargin < 2
  radius = 1;
end

if nargin < 3
  center = [];
end
   
if nargin < 5
  sang = 0;
  eang = 2*pi;
end

surf = nrlextrude(nrlcirc(radius,center,sang,eang),[0.0 0.0 height]);

end

%% demo
%  srf = nrlcylind(3,1,[],deg2rad(270),deg2rad(180));
%  nrlplot(srf,[20,20]);
%  axis equal;
%  title('Cylinderical section by extrusion of a circular arc.');
%  hold off



