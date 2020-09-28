function nrblplot (nurl, subd, varargin)
% 
% NRBLPLOT: Plot a NURL curve or surface, or the boundary of a NURL volume.
% 
% Calling Sequence:
% 
%   nrblplot (nurl, subd)
%   nrblplot (nurl, subd, p, v)
% 
% INPUT:
% 
%   nurl		: NURL curve, surface or volume, see nrbmak.
% 
%   npnts	: Number of evaluation points, for a surface or volume, a row 
%       vector with the number of points along each direction.
% 
%   [p,v]       : property/value options
%
%               Valid property/value pairs include:
%
%               Property        Value/{Default}
%               -----------------------------------
%               light           {off} | on
%               colormap        {'copper'}
%
% Example:
%
%   Plot the test surface with 20 points along the U direction
%   and 30 along the V direction
%
%   nrblplot(nrltestsrf, [20 30])
%

nargs = nargin;
if nargs < 2
  error ('Need a NURL to plot and the number of subdivisions!');
elseif rem(nargs+2,2)
  error ('Param value pairs expected')
end

% Default values
light='off';
cmap='summer';

% Recover Param/Value pairs from argument list
for i=1:2:nargs-2
  Param = varargin{i};
  Value = varargin{i+1};
  if (~ischar (Param))
    error ('Parameter must be a string')
  elseif size(Param,1)~=1
    error ('Parameter must be a non-empty single row string.')
  end
  switch lower (Param)
  case 'light'
    light = lower (Value);
    if (~ischar (light))
      error ('light must be a string.')
    elseif ~(strcmp(light,'off') || strcmp(light,'on'))
      error ('light must be off | on')
    end
  case 'colormap'
    if ischar (Value)
      cmap = lower(Value);
    elseif size (Value, 2) ~= 3
      error ('colormap must be a string or have exactly three columns.')
    else
      cmap=Value;
    end
  otherwise
    error ('Unknown parameter: %s', Param)
  end
end

colormap (cmap);

% convert the number of subdivisions in number of points
subd = subd+1;

% plot the curve or surface
if (iscell (nurl.knots))
 if (size (nurl.knots,2) == 2) % plot a NURL surface
  p = nrleval (nurl, {linspace(0, 1, subd(1)) ...
                       linspace(0, 1, subd(2))});
  if (strcmp (light,'on'))
    % light surface
    surfl (squeeze(p(1,:,:)), squeeze(p(2,:,:)), squeeze(p(3,:,:)));
    shading interp;
  else 
    surf (squeeze (p(1,:,:)), squeeze (p(2,:,:)), squeeze (p(3,:,:)));
    shading faceted;
  end
 elseif (size (nurl.knots,2) == 3) % plot the boundaries of a NURL volume
  bnd = nrlextract (nurl);
  hold_flag = ishold;
  nrblplot (bnd(1), subd(2:3), varargin{:});
  hold on
  nrblplot (bnd(2), subd(2:3), varargin{:});
  nrblplot (bnd(3), subd([1 3]), varargin{:});
  nrblplot (bnd(4), subd([1 3]), varargin{:});
  nrblplot (bnd(5), subd(1:2), varargin{:});
  nrblplot (bnd(6), subd(1:2), varargin{:});
  
  if (~hold_flag)
    hold off
  end
 
 else
  error ('nrbplot: some argument is not correct')
 end
else
  % plot a NURL curve
  p = nrleval (nurl, linspace (0, 1, subd));

  if (any (nurl.coefs(3,:)))
    % 3D curve
    plot3 (p(1,:), p(2,:), p(3,:)); 
    grid on;
  else
    % 2D curve
    plot (p(1,:), p(2,:));
  end
end
axis equal;

end

%% demo - curve
%  crv = nrltestcrv;
%  nrblplot(crv,100)
%  title('Test curve')
%  hold off

%% demo 
%  coefs = [0.0 7.5 15.0 25.0 35.0 30.0 27.5 30.0;
%           0.0 2.5  0.0 -5.0  5.0 15.0 22.5 30.0];
%  knots = [0.0 0.0 0.0 1/6 1/3 1/2 2/3 5/6 1.0 1.0 1.0];
% 
%  geom = [
%  nrb2nrl(nrbmak(coefs,knots))
%  nrlline([30.0 30.0],[20.0 30.0])
%  nrlline([20.0 30.0],[20.0 20.0])
%  nrlcirc(10.0,[10.0 20.0],1.5*pi,0.0)
%  nrlline([10.0 10.0],[0.0 10.0])
%  nrlline([0.0 10.0],[0.0 0.0])
%  nrlcirc(5.0,[22.5 7.5])
%  ];
% 
%  ng = length(geom);
%  for i = 1:ng
%    nrblplot(geom(i),500);
%    hold on;
%  end
%  hold off;
%  axis equal;
%  title('2D Geometry formed by a series of NURL curves');

%% demo
%  sphere = nrlrevolve(nrlcirc(1,[],0.0,pi),[0.0 0.0 0.0],[1.0 0.0 0.0]);
%  nrblplot(sphere,[40 40],'light','on');
%  title('Ball and torus - surface construction by revolution');
%  hold on;
%  torus = nrlrevolve(nrlcirc(0.2,[0.9 1.0]),[0.0 0.0 0.0],[1.0 0.0 0.0]);
%  nrblplot(torus,[40 40],'light','on');
%  hold off
