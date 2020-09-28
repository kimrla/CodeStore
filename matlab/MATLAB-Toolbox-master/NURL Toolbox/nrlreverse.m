function nrl = nrlreverse(nrl, idir)
%
% NRLREVERSE: Reverse the evaluation directions of a NURL geometry.
% 
% Calling Sequence:
% 
%   rnrl = nrlreverse(nrl);
%   rnrl = nrlreverse(nrl, idir);
% 
% INPUT:
% 
%   nrl		: NURL data structure, see nrbmak.
%   idir        : vector of directions to reverse.
%
% OUTPUT:
% 
%   rnrl	: Reversed NURL.
% 
% Description:
% 
%   Utility function to reverse the evaluation direction of a NURL
%   curve or surface.
%

if (nargin > 2)
  error('Incorrect number of input arguments');
end

if (iscell(nrl.knots))
  % reverse a NURL surface or volume
  ndim = numel (nrl.knots);
  if (nargin == 1 || isempty (idir))
    idir = 1:ndim;
  end
  for ii = idir
    nrl.knots{ii} = sort (nrl.knots{ii}(end) - nrl.knots{ii});
    nrl.intervals{ii} = sort (nrl.intervals{ii}(end) - nrl.intervals{ii});
    nrl.coefs = flip (nrl.coefs, ii+1);
  end

else
  % reverse a NURL curve
  nrl.knots = sort (nrl.knots(end) - nrl.knots);
  nrl.intervals = sort (nrl.intervals(end) - nrl.intervals); 
  nrl.coefs = fliplr (nrl.coefs);
end

end

%% demo - curve
%  pnts = [0.5 1.5 3.0 7.5 8.5;
%          3.0 5.5 1.5 4.0 4.5;
%          0.0 0.0 0.0 0.0 0.0];
%  crv1 = nrbmak(pnts,[0 0 0 1/2 3/4 1 1 1]);
%  crv1=nrb2nrl(crv1);
%  crv2 = nrlreverse(crv1);
%  nrblplot(crv1,100)
%  hold on
%  nrblplot(crv2,100)
%  title('The curve and its reverse are the same')
%  hold off

%% demo - surface
%  srf  = nrlrevolve(nrlline([1 0],[2 0]), [0 0 0], [0 0 1], pi/2);
%  srf  = nrlintins (srf, {0.3, 0.6});
%  srf2 = nrlreverse (srf);
%  nrlctrlplot(srf2);



