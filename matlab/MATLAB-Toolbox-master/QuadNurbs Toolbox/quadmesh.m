function hh = quadmesh(quad, varargin)
%QUADMESH Triangular mesh plot
%   QUADMESH(QUAD,X,Y,Z,C) displays the quadrangles defined in the M-by-3
%   face matrix TRI as a mesh.  A row of QUAD contains indexes into
%   the X,Y, and Z vertex vectors to define a single quadrangular face.
%   The edge color is defined by the vector C.
%
%   QUADMESH(QUAD,X,Y,Z) uses C = Z, so color is proportional to surface
%   height.
%
%   QUADMESH(QUAD,X,Y) displays the quadrangles in a 2-d plot.
%
%   H = QUADMESH(...) returns a handle to the displayed quadrangles.
%
%   QUADMESH(...,'param','value','param','value'...) allows additional
%   patch param/value pairs to be used when creating the patch object. 
%

ax = axescheck(varargin{:});
ax = newplot(ax);

if nargin == 3 || (nargin > 4 && ischar(varargin{3}))
  d = quad(:,[1 2 3 4 1])';
  x = varargin{1};
  y = varargin{2};
  if nargin == 3
    h = plot(ax, x(d), y(d));
  else
    z = varargin{3};  
    h = plot(ax, x(d), y(d),z,varargin{4},varargin{5:end});
  end
  if nargout == 1, hh = h; end
  return;
else
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    trids = quad;
    if nargin>4 && rem(nargin-4,2)==1
      c = varargin{4};
      start = 5;
    elseif nargin<3
      error(message('MATLAB:trimesh:NotEnoughInputs'));
    else
      start = 4;
      c = z;
    end
end

if ischar(get(ax,'color')),
  fc = get(gcf,'Color');
else
  fc = get(ax,'color');
end

h = patch('faces',trids,'vertices',[x(:) y(:) z(:)],'facevertexcdata',c(:),...
	  'facecolor',fc,'edgecolor',get(ax,'DefaultSurfaceFaceColor'),...
	  'facelighting', 'none', 'edgelighting', 'flat',...
      'parent',ax,...
	  varargin{start:end});
if ~ishold(ax), view(ax,3), grid(ax,'on'), end
if nargout == 1, hh = h; end









