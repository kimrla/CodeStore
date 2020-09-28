function nurls = nrltform(nurls,tmat)
% 
% NRBTFORM: Apply transformation matrix to the NURLS.
% 
% Calling Sequence:
% 
%   tnurls = nrbtform(nurls,tmatrix);
% 
% INPUT:
% 
%   nurbs	: NURLS data structure (see nrbmak for details).
% 
%   tmatrix     : Transformation matrix, a matrix of size (4,4) defining
%                 a single or multiple transformations.
%
% OUTPUT:
%
%   tnurbs	: The return transformed NURLS data structure.
% 
% Description:
% 
%   The NURLS is transform as defined a transformation matrix of size (4,4),
%   such as a rotation, translation or change in scale. The transformation
%   matrix can define a single transformation or multiple series of
%   transformations. The matrix can be simply constructed by the functions
%   vecscale, vectrans and vecrot, and also vecrotx, vecroty, and vecrotz.
%     
% Examples:
% 
%   Rotate a square by 45 degrees about the z axis.
%
%   rsqr = nrltform(nrlrect(), vecrotz(deg2rad(45)));
%   nrlplot(rsqr, 1000);
% 
% See also:
% 
%   vecscale, vectrans, vecrot, vecrotx, vecroty, vecrotz
%

if nargin < 2
  error('Not enough input arguments!');
end;

nurls.coefs(1:3, :)=nurls.coefs(1:3, :).*repmat(nurls.coefs(4, :), 3, 1);

if iscell(nurls.knots)
 if size(nurls.knots,2) == 2
  % NURBS is a surface
  [dim,nu,nv] = size(nurls.coefs);
  nurls.coefs = reshape(tmat*reshape(nurls.coefs,dim,nu*nv),[dim nu nv]);
 elseif size(nurls.knots,2) == 3
  % NURBS is a volume
  [dim,nu,nv,nw] = size(nurls.coefs);
  nurls.coefs = reshape(tmat*reshape(nurls.coefs,dim,nu*nv*nw),[dim nu nv nw]);
 end
else
  % NURBS is a curve
  nurls.coefs = tmat*nurls.coefs;
end

nurls.coefs(1:3, :)=nurls.coefs(1:3, :)./repmat(nurls.coefs(4, :), 3, 1);
end

%!demo
%! xx = vectrans([2.0 1.0])*vecroty(pi/8)*vecrotx(pi/4)*vecscale([1.0 2.0]);
%! c0 = nrltform(nrlcirc, xx);
%! nrlplot(c0,50);
%! grid on
%! title('Construction of an ellipse by transforming a unit circle.');
%! hold off



