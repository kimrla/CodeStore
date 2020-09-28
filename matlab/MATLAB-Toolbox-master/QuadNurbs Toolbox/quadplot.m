function hh = quadplot(qnrb,varargin)
%QUADPLOT Plots a 2D quadrangulation
%   QUADPLOT(QNRB) displays the quadrangles defined in the
%   QNRB parametric plane.  
%
%   QUADPLOT(QUAD, COLOR) uses the string COLOR as the line color.
%
%   H = QUADPLOT(...) returns a line handle representing the displayed
%   quadrangles edges.
%
%   QUADPLOT(...,'param','value','param','value'...) allows additional
%   line param/value pairs to be used when creating the plot.
%

narginchk(1,inf);

x = qnrb.nodes(:,1);
y = qnrb.nodes(:,2);
edges = qnrb.qedges;
if (nargin == 1) || (mod(nargin-1,2) == 0)
  c = 'blue';
  start = 1;
else 
  c = varargin{1};
  start = 2;
end

x = x(edges)';
y = y(edges)';
nedges = size(x,2);
x = [x; NaN(1,nedges)];
y = [y; NaN(1,nedges)];
x = x(:);
y = y(:);
h = plot(x,y,c,varargin{start:end});
if nargout == 1, hh = h; end











