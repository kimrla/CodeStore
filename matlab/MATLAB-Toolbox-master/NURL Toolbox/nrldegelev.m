function inurls = nrldegelev(nurls, ntimes)
% 
% NRLDEGELEV: Elevate the degree of the NURL curve, surface or volume.
% 
% Calling Sequence:
% 
%   ecrv = nrldegelev(crv, utimes);
%   esrf = nrldegelev(srf, [utimes,vtimes]);
%   evol = nrldegelev(vol, [utimes, vtimes, wtimes]);
% 
% INPUT:
% 
%   crv		: NURL curve, see nrbmak.
% 
%   srf		: NURL surface, see nrbmak.
% 
%   vol		: NURL volume, see nrbmak.
% 
%   utimes	: Increase the degree along U direction utimes.
% 
%   vtimes	: Increase the degree along V direction vtimes.
% 
%   wtimes	: Increase the degree along W direction vtimes.
%
% OUTPUT:
%
%   ecrv	: new NURL structure for a curve with degree elevated.
% 
%   esrf	: new NURL structure for a surface with degree elevated.
% 
%   evol	: new NURL structure for a volume with degree elevated.
% 
% 
% Description:
% 
%   Degree elevates the NURL curve or surface. This function uses the
%   NUL function nuldegelev.
% 
% Examples:
% 
%   Increase the NURL surface twice along the V direction.
%   esrf = nrldegelev(srf, [0, 2]); 
% 
% See also:
% 
%   nuldegelev
%

if nargin < 2
  error('Input argument must include the NURL and degree increment.');
end

if ~isstruct(nurls)
  error('NURL representation is not structure!');
end

if ~strcmp(nurls.form,'L-NURL')
  error('Not a recognised NURL representation');
end

degree = nurls.order;
p=ntimes<0; 
ntimes(p)=0; 
ntimes=round(ntimes);
if length(degree)~=length(ntimes)
    error('Not correct input of degrees to be elevated.');
end

% Transform NURL coefs into homogeneous coordinates (wx,wy,wz)
nurls.coefs(1:3, :)=nurls.coefs(1:3, :).*repmat(nurls.coefs(4, :), 3, 1);

% Degree elevation
if iscell(nurls.knots)
 if size(nurls.knots,2) == 3
  % NURBS represents a volume
  [dim,num1,num2,num3] = size(nurls.coefs);

  % Degree elevate along the w direction
  if ntimes(3) == 0
    coefs = nurls.coefs;
    knots{3} = nurls.knots{3};
  else
    coefs = reshape(nurls.coefs,4*num1*num2,num3);
    [coefs,knots{3}] = nuldegelev(degree(3),coefs,nurls.knots{3},nurls.intervals{3},ntimes(3));
    num3 = size(coefs,2);
    coefs = reshape(coefs,[4 num1 num2 num3]);
  end

  % Degree elevate along the v direction
  if ntimes(2) == 0
    knots{2} = nurls.knots{2};
  else
    coefs = permute(coefs,[1 2 4 3]);
    coefs = reshape(coefs,4*num1*num3,num2);
    [coefs,knots{2}] = nuldegelev(degree(2),coefs,nurls.knots{2},nurls.intervals{2},ntimes(2));
    num2 = size(coefs,2);
    coefs = reshape(coefs,[4 num1 num3 num2]);
    coefs = permute(coefs,[1 2 4 3]);
  end

  % Degree elevate along the u direction
  if ntimes(1) == 0
    knots{1} = nurls.knots{1};
  else
    coefs = permute(coefs,[1 3 4 2]);
    coefs = reshape(coefs,4*num2*num3,num1);
    [coefs,knots{1}] = nuldegelev(degree(1),coefs,nurls.knots{1},nurls.intervals{1},ntimes(1));
    coefs = reshape(coefs,[4 num2 num3 size(coefs,2)]); 
    coefs = permute(coefs,[1 4 2 3]); 
  end 

 elseif size(nurls.knots,2) == 2
  % NURL represents a surface
  [dim,num1,num2] = size(nurls.coefs);

  % Degree elevate along the v direction
  if ntimes(2) == 0
    coefs = nurls.coefs;
    knots{2} = nurls.knots{2};
  else
    coefs = reshape(nurls.coefs,4*num1,num2);
    [coefs,knots{2}] = nuldegelev(degree(2),coefs,nurls.knots{2},nurls.intervals{2},ntimes(2));
    num2 = size(coefs,2);
    coefs = reshape(coefs,[4 num1 num2]);
  end

  % Degree elevate along the u direction
  if ntimes(1) == 0
    knots{1} = nurls.knots{1};
  else
    coefs = permute(coefs,[1 3 2]);
    coefs = reshape(coefs,4*num2,num1);
    [coefs,knots{1}] = nuldegelev(degree(1),coefs,nurls.knots{1},nurls.intervals{1},ntimes(1));
    coefs = reshape(coefs,[4 num2 size(coefs,2)]);
    coefs = permute(coefs,[1 3 2]);
  end 
 end
else

  % NURL represents a curve
  if isempty(ntimes) || ntimes==0
    coefs = nurls.coefs;
    knots = nurls.knots;
  else
    [coefs,knots] = nuldegelev(degree,nurls.coefs,nurls.knots,nurls.intervals,ntimes);
  end
  
end

% Construct new NURLS
coefs(1:3, :)=coefs(1:3, :)./repmat(coefs(4, :), 3, 1);
inurls = nrlmake(coefs,knots,nurls.intervals,degree+ntimes);

end

%% demo
%  crv = nrbtestcrv;
%  crv=nrb2nrl(crv);
%  plot(crv.coefs(1,:),crv.coefs(2,:),'bo')
%  title('Degree elevation along test curve: curve and control points.');
%  hold on;
%  nrlplot(crv,49);
%  
%  icrv = nrldegelev(crv, 1);
%  
%  plot(icrv.coefs(1,:),icrv.coefs(2,:),'ro')
%  
%  hold off;
