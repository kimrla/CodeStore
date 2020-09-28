function volcplot(x, y, z, c, varargin)

% VOLCPLOT: Plot a volume with colored data
% 
% Calling Sequence:
% 
%   volplot(x, y, z)
% 
%   volplot(x, y, z, c)
% 
%   volplot(x, y, z, c, [p,v])
% 
% INPUT:
% 
%   x, y, z    : 3D coordinates
%
%    c    :  values on (x, y, z), uses c to define color
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

nargs = nargin;
if nargs == 3
    volplot(x, y, z);
    return;
elseif nargs < 3
  error ('Need a volume (x, y, z) or (x, y, z, c) to plot!');
elseif rem(nargs+4,2)
  error ('Param value pairs expected')
end

% Default values
light='off';
cmap='summer';

% Recover Param/Value pairs from argument list
for i=1:2:nargs-4
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
    elseif ~(strcmp(light,'off') | strcmp(light,'on'))
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

hold_flag = ishold;
if (strcmp (light,'on'))
    surf(squeeze(x(1,:,:)), squeeze(y(1,:,:)), squeeze(z(1,:,:)), squeeze(c(1,:,:)));
    hold on
    surf(squeeze(x(end,:,:)), squeeze(y(end,:,:)), squeeze(z(end,:,:)), squeeze(c(end,:,:)));
    surf(squeeze(x(:,1,:)), squeeze(y(:,1,:)), squeeze(z(:,1,:)), squeeze(c(:,1,:)));
    surf(squeeze(x(:,end,:)), squeeze(y(:,end,:)), squeeze(z(:,end,:)), squeeze(c(:,end,:)));
    surf(squeeze(x(:,:,1)), squeeze(y(:,:,1)), squeeze(z(:,:,1)), squeeze(c(:,:,1)));
    surf(squeeze(x(:,:,end)), squeeze(y(:,:,end)), squeeze(z(:,:,end)), squeeze(c(:,:,end)));
    shading interp;
else
    surf(squeeze(x(1,:,:)), squeeze(y(1,:,:)), squeeze(z(1,:,:)), squeeze(c(1,:,:)));
    hold on
    surf(squeeze(x(end,:,:)), squeeze(y(end,:,:)), squeeze(z(end,:,:)), squeeze(c(end,:,:)));
    surf(squeeze(x(:,1,:)), squeeze(y(:,1,:)), squeeze(z(:,1,:)), squeeze(c(:,1,:)));
    surf(squeeze(x(:,end,:)), squeeze(y(:,end,:)), squeeze(z(:,end,:)), squeeze(c(:,end,:)));
    surf(squeeze(x(:,:,1)), squeeze(y(:,:,1)), squeeze(z(:,:,1)), squeeze(c(:,:,1)));
    surf(squeeze(x(:,:,end)), squeeze(y(:,:,end)), squeeze(z(:,:,end)), squeeze(c(:,:,end)));
    shading faceted;
end

if (~hold_flag)
    hold off
end


