function srf = nrlruled (crv1, crv2)

% NRLRULED: Construct a ruled surface between two NURL curves.
% 
% Calling Sequence:
% 
%   srf = nrlruled(crv1, crv2)
% 
% INPUT:
% 
%   crv1	: First NURL curve, see nrbmak.
% 
%   crv2	: Second NURL curve, see nrbmak.
%
% OUTPUT:
% 
%   srf		: Ruled NURL surface.
% 
% Description:
% 
%   Constructs a ruled surface between two NURL curves. The ruled surface is
%   ruled along the V direction.
% 
% Examples:
% 
%   Construct a ruled surface between a semicircle and a straight line.
% 
%   cir = nrbcirc(1,[0 0 0],0,pi);
%   line = nrbline([-1 0.5 1],[1 0.5 1]);
%   srf = nrbruled(cir,line);
%   nrbplot(srf,[20 20]);
%

if (iscell(crv1.knots) || iscell(crv2.knots))
  error ('Both NURL must be curves');
end

% ensure both curves have a common degree
d = max ([crv1.order, crv2.order]);
crv1 = nrldegelev (crv1, d - crv1.order);
crv2 = nrldegelev (crv2, d - crv2.order);

% merge the knot vectors, to obtain a common knot vector
k1 = crv1.intervals;
k2 = crv2.intervals;
ku = unique ([k1 k2]);
crv1 = nrlintins (crv1, ku);
crv2 = nrlintins (crv2, ku);
nk1=countknts(crv1.intervals, crv1.knots); 
nk2=countknts(crv2.intervals, crv2.knots); 
nk=max([nk1, nk2]);
in1=nk-nk1;
in2=nk-nk2;
crv1 = nrlkntins (crv1, in1);
crv2 = nrlkntins (crv2, in2);

coefs(:,:,1) = crv1.coefs;
coefs(:,:,2) = crv2.coefs;
srf = nrlmake (coefs, {crv1.knots, [0 1]}, {crv1.intervals, [0 1]}, [crv1.order, 1]);

end


%% demo
% pnts = [0.5 1.5 4.5 3.0 7.5 6.0 8.5;
%         3.0 5.5 5.5 1.5 1.5 4.0 4.5;
%         0.0 0.0 0.0 0.0 0.0 0.0 0.0];
% crv1 = nrbmak (pnts,[0 0 0 1/4 1/2 3/4 3/4 1 1 1]);
% crv2 = nrbtform (nrbcirc (4,[4.5;0],pi,0.0),vectrans([0.0 4.0 -4.0]));
% crv1=nrb2nrl(crv1); crv2=nrb2nrl(crv2);
% srf = nrlruled (crv1,crv2);
% nrlplot (srf,[41 20]);
% title ('Ruled surface construction from two NURL curves.');
% hold off




