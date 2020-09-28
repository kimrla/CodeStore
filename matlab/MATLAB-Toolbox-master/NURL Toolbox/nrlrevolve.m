function surf = nrlrevolve(curve, pnt, vec, theta)

% 
% NRLREVOLVE: Construct a NURL surface by revolving a NURL curve, or
%  construct a NURL volume by revolving a NURL surface.
% 
% Calling Sequence:
% 
%   srf = nrlrevolve(crv,pnt,vec,ang)
% 
% INPUT:
% 
%   crv		: NURL curve or surface to revolve, see nrbmak.
% 
%   pnt		: Coordinates of the point used to define the axis
%               of rotation.
% 
%   vec		: Vector defining the direction of the rotation axis.
% 
%   ang		: Angle to revolve the curve, default 2*pi
%
% OUTPUT:
%
%   srf		: constructed surface or volume
% 
% Description:
% 
%   Construct a NURL surface by revolving the profile NURL curve around
%   an axis defined by a point and vector.
% 
% Examples:
% 
%   Construct a sphere by rotating a semicircle around a x-axis.
%
%   crv = nrlcirc(1.0,[0 0 0],0,pi);
%   srf = nrlrevolve(crv,[0 0 0],[1 0 0]);
%   nrlplot(srf,[20 20]);
%
% NOTE:
%
%   The algorithm:
%
%     1) vectrans the point to the origin (0,0,0)
%     2) rotate the vector into alignment with the z-axis
%         for each control point along the curve
%     3) determine the radius and angle of control
%         point to the z-axis
%     4) construct a circular arc in the x-y plane with 
%         this radius and start angle and sweep angle theta 
%     5) combine the arc and profile, coefs and weights.
%         next control point
%     6) rotate and vectrans the surface back into position
%          by reversing 1 and 2.
%

if (nargin < 3)
  error('Not enough arguments to construct revolved surface');
end

if (nargin < 4)
  theta = 2.0*pi;
end

if (iscell (curve.knots) && numel(curve.knots) == 3)
  error('The function nrbrevolve is not yet ready to create volumes') 
end

% Translate curve the center point to the origin
if isempty(pnt)
  pnt = zeros(3,1);
end

if length(pnt) ~= 3
  error('All point and vector coordinates must be 3D');
end

% Translate and rotate the original curve or surface into alignment with the z-axis
T  = vectrans(-pnt);
angx = vecangle(vec(1),vec(3));
RY = vecroty(-angx);
vectmp = RY*[vecnorm(vec(:));1.0];
angy = vecangle(vectmp(2),vectmp(3));
RX = vecrotx(angy);
curve = nrltform(curve,RX*RY*T);

% Construct an arc 
arc = nrlcirc(1.0,[],0.0,theta);
curve.coefs(1:3, :)=curve.coefs(1:3, :).*repmat(curve.coefs(4, :), 3, 1);
arc.coefs(1:3, :)=arc.coefs(1:3, :).*repmat(arc.coefs(4, :), 3, 1);

if (iscell (curve.knots))
% Construct the revolved volume
  coefs = zeros([4 arc.number curve.number]);
  angle = squeeze (vecangle(curve.coefs(2,:,:),curve.coefs(1,:,:)));
  radius = squeeze (vecmag(curve.coefs(1:2,:,:)));
  for i = 1:curve.number(1)  
    for j = 1:curve.number(2)  
      coefs(:,:,i,j) = vecrotz(angle(i,j))*vectrans([0.0 0.0 curve.coefs(3,i,j)])*...
          vecscale([radius(i,j) radius(i,j)])*arc.coefs;
      coefs(4,:,i,j) = coefs(4,:,i,j)*curve.coefs(4,i,j);
   end
  end
  surf = nrlmake(coefs,{arc.knots, curve.knots{:}},{arc.intervals, curve.intervals{:}},[arc.order, curve.order]);
else
% Construct the revolved surface
  coefs = zeros(4, arc.number, curve.number);
  angle = vecangle(curve.coefs(2,:),curve.coefs(1,:));
  radius = vecmag(curve.coefs(1:2,:));
  for i = 1:curve.number  
    coefs(:,:,i) = vecrotz(angle(i))*vectrans([0.0 0.0 curve.coefs(3,i)])*...
          vecscale([radius(i) radius(i)])*arc.coefs;
    coefs(4,:,i) = coefs(4,:,i)*curve.coefs(4,i);
  end
  surf = nrlmake(coefs,{arc.knots, curve.knots},{arc.intervals, curve.intervals},[arc.order, curve.order]);
end
surf.coefs(1:3, :)=surf.coefs(1:3, :)./repmat(surf.coefs(4, :), 3, 1);

% Rotate and vectrans the surface back into position
T = vectrans(pnt);
RX = vecrotx(-angy);
RY = vecroty(angx);
surf = nrltform(surf,T*RY*RX);  

end

%% demo
% sphere = nrlrevolve(nrlcirc(1,[],0.0,pi),[0.0 0.0 0.0],[1.0 0.0 0.0]);
% nrlplot(sphere,[40 40]);
% title('Ball and tori - surface construction by revolution');
% hold on;
% torus = nrlrevolve(nrlcirc(0.2,[0.9 1.0]),[0.0 0.0 0.0],[1.0 0.0 0.0]);
% nrlplot(torus,[40 20]);
% nrlplot(nrltform(torus,vectrans([-1.8])),[40 20]);
% hold off;

%% demo
% pnts = [3.0 5.5 5.5 1.5 1.5 4.0 4.5;
%         0.0 0.0 0.0 0.0 0.0 0.0 0.0;
%         0.5 1.5 4.5 3.0 7.5 6.0 8.5];
% crv = nrbmak(pnts,[0 0 0 1/4 1/2 3/4 3/4 1 1 1]);
% crv=nrb2nrl(crv);
% 
% xx = vecrotz(deg2rad(25))*vecroty(deg2rad(15))*vecrotx(deg2rad(20));
% nrl = nrltform(crv,vectrans([5 5])*xx);
% 
% pnt = [5 5 0]';
% vec = xx*[0 0 1 1]';
% srf = nrlrevolve(nrl,pnt,vec(1:3));
% 
% p = nrleval(srf,{linspace(0.0,1.0,100) linspace(0.0,1.0,100)});
% surfl(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)));
% title('Construct of a 3D surface by revolution of a curve.');
% shading interp;
% colormap(copper);
% axis equal;
% hold off

%% demo
%  crv1 = nrlcirc(1,[0 0],0, pi/2);
%  crv2 = nrlcirc(2,[0 0],0, pi/2);
%  srf = nrlruled (crv1, crv2);
%  srf = nrltform (srf, [1 0 0 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]);
%  vol = nrlrevolve (srf, [0 0 0], [1 0 0], pi/2);
%  nrblplot(vol, [30 30 30], 'light', 'on')







