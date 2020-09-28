function tvol = nrlpermute (vol, ord)
% 
% NRLPERMUTE: Rearrange the directions of a NURL volume or surface.
% 
% Calling Sequence:
% 
%   tvol = nrlpermute(vol,order)
%
% INPUT:
% 
%   vol		: NURL volume or surface, see nrlmake.
%   order   : the order to rearrange the directions of the NURL entity.
%
% OUTPUT:
% 
%   tvol	: NURL volume or surface with rearranged directions.
% 
% Description:
% 
%   Utility function that rearranges the directions of a NURL volume or
%   surface. For surfaces, nrlpermute(srf,[1 2]) is the same as
%   nrltransp(srf). NURL curves cannot be rearranged.
%
% Example:
%
%    nrlpermute (vol, [1 3 2])
%

if (~iscell(vol.knots))
  error('A NURL curve cannot be rearranged.');
end

tvol = nrlmake (permute (vol.coefs, [1, ord+1]), {vol.knots{ord}}, {vol.intervals{ord}}, vol.order(ord));

%% demo
%! vol = nrlrevolve (nrl4surf ([1 0], [2 0], [1 1], [2 1]), [0 0 0], [0 1 0], pi/8);
%! nrblplot(vol,[5 10 20]);
%! title('NURL volume and the same after reordering the directions')
%! hold on
%! vol.coefs(1,:,:) = vol.coefs(1,:,:) + 2;
%! vol = nrlpermute(vol,[2 3 1]);
%! nrblplot(vol,[5 10 20]);
%! hold off

