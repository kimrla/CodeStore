function hh =quadsurf(quad, varargin)
%QUADSURF Triangular surface plot
%   QUADSURF(QUAD,X,Y,Z,C) displays the quadrangles defined in the M-by-3
%   face matrix QUAD as a surface.  A row of QUAD contains indexes into
%   the X,Y, and Z vertex vectors to define a single quadrangular face.
%   The color is defined by the vector C.
%
%   QUADSURF(QUAD,X,Y,Z) uses C = Z, so color is proportional to surface
%   height.
%
%   H = QUADSURF(...) returns a patch handle.
%
%   QUADSURF(...,'param','value','param','value'...) allows additional
%   patch param/value pairs to be used when creating the patch object. 
%

narginchk(1,inf);

ax = axescheck(varargin{:});
ax = newplot(ax);

x = varargin{1};
y = varargin{2};
z = varargin{3};
trids = quad;
if nargin>4 && rem(nargin-4,2)==1, 
    c = varargin{4};
    start = 5;
else
    c = z;
    start = 4;
end
h = patch('faces',trids,'vertices',[x(:) y(:) z(:)],'facevertexcdata', c(:),...
      'facecolor',get(ax,'DefaultSurfaceFaceColor'), ...
      'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),'parent',ax,...
      varargin{start:end});
if ~ishold(ax), view(ax,3), grid(ax,'on'), end
if nargout==1, hh = h; end








