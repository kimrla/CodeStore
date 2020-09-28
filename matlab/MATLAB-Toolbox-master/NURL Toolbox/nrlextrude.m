function srf = nrlextrude(curve,vector)

% NRLEXTRUDE: Construct a NURL surface by extruding a NURL curve, or 
%  construct a NURL volume by extruding a NURL surface.
% 
% Calling Sequence:
% 
%   srf = nrbextrude(crv,vec);
% 
% INPUT:
% 
%   crv		: NURL curve or surface to extrude, see nrbmak.
% 
%   vec		: Vector along which the entity is extruded.
%
% OUTPUT: 
% 
%   srf		: NURL surface or volume constructed.
% 
% Description:
% 
%   Constructs either a NURL surface by extruding a NURL curve along a  
%   defined vector, or a NURL volume by extruding a NURL surface. In the 
%   first case, the NURL curve forms the U direction of the surface edge, and
%   is extruded along the vector in the V direction. In the second case, the 
%   original surface forms the U and V direction of the volume, and is extruded
%   along the W direction.
%
% Examples:
% 
%   Form a hollow cylinder by extruding a circle along the z-axis.
%
%   srf = nrbextrude(nrbcirc, [0,0,1]);
%


if (nargin < 2)
  error('Error too few input arguments!');
end

% Transform NURL coefs into homogeneous coordinates (wx,wy,wz)
curve.coefs(1:3, :)=curve.coefs(1:3, :).*repmat(curve.coefs(4, :), 3, 1);

if (iscell (curve.knots))
      if (numel (curve.knots) == 3)
        error('Nurl volumes cannot be extruded!');
      end
      for ii = 1:size(curve.coefs,3)
        coefs(:,:,ii) = vectrans(vector) * squeeze (curve.coefs(:,:,ii));
      end
      coefs = cat(4,curve.coefs,coefs);
      srf = nrlmake(coefs, {curve.knots{:}, [0 1]}, {curve.intervals{:}, [0 1]}, [curve.order, 1]);
else
      coefs = cat(3,curve.coefs,vectrans(vector)*curve.coefs);
      srf = nrlmake(coefs,{curve.knots, [0 1]}, {curve.intervals, [0 1]}, [curve.order, 1]);
end
srf.coefs(1:3, :)=srf.coefs(1:3, :)./repmat(srf.coefs(4, :), 3, 1);
end

%% demo - curve
% crv = nrbtestcrv;
% crv=nrb2nrl(crv);
% srf = nrlextrude(crv,[0 0 5]);
% nrlplot(srf,[40 10]);
% title('Extrusion of a test curve along the z-axis');
% hold off
%
%% demo - surface
% crv1 = nrlcirc (1, [0 0], 0, pi/2);
% crv2 = nrlcirc (2, [0 0], 0, pi/2);
% srf  = nrlruled (crv1, crv2);
% vol  = nrlextrude (srf, [0 0 1]);
% nrlplot (vol, [30 10 10])
% title ('Extrusion of the quarter of a ring')
%
%% demo - surface
% srf = nrltestsrf;
% vol = nrlextrude(srf, [0 0 10]);
% nrlplot(vol,[20 20 20]);
% title('Extrusion of a test surface along the z-axis');
% hold off
