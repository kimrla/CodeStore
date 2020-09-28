function inurl = nrlintins(nurls, intins)
% 
% NRLINTINS: Insert a single or multiple intervals into a NURL curve,
%            surface or volume.
% 
% Calling Sequence:
% 
%   icrv = nrlkntins(crv,iuintvs);
%   isrf = nrlkntins(srf,{iuintvs ivintvs});
%   ivol = nrlkntins(vol,{iuintvs ivintvs iwintvs});
% 
% INPUT:
% 
%   crv		: NURL curve, see nrbmak.
% 
%   srf		: NURL surface, see nrbmak.
% 
%   srf		: NURL volume, see nrbmak.
% 
%   iuintvs	: Intervals to be inserted along U direction.
% 
%   ivintvs	: Intervals to be inserted along V direction.
% 
%   iwintvs	: Intervals to be inserted along W direction.
% 
% OUTPUT:
% 
%   icrv	: new NURL structure for a curve with intervals inserted.
% 
%   isrf	: new NURL structure for a surface with intervals inserted.
% 
%   ivol	: new NURL structure for a volume with intervals inserted.
% 
% Description:
% 
%   Inserts intervals into the NURL data structure. Intervals along 
%   the V direction can only inserted into NURL surfaces, not 
%   curve that are always defined along the U direction. 
%   This function use the NUL function nulintins. 
%   Notice: Insertions outside the interval of definition will be ignored
% 
% Examples:
% 
%   Insert two intervals into a curve, one at 0.3 and another at 0.4
% 
%   icrv = nrlintins(crv, [0.3 0.4])
% 
%   Insert into a surface two intervals as (1) into the U interval
%   sequence and one interval into the V interval sequence at 0.5.
%
%   isrf = nrlintins(srf, {[0.3 0.4] [0.5]})
% 
% See also:
% 
%   nulintins
%


if nargin < 2
  error('Input argument must include the NURL and knots to be inserted');
end

if ~isstruct(nurls)
  error('NURL representation is not structure!');
end

if ~strcmp(nurls.form,'L-NURL')
  error('Not a recognised NURL representation');
end

degree = nurls.order;
nurls.coefs(1:3, :)=nurls.coefs(1:3, :).*repmat(nurls.coefs(4, :), 3, 1);

if iscell(nurls.knots)
     if size(nurls.knots,2)==3
      % NURL represents a volume
      num1 = nurls.number(1);
      num2 = nurls.number(2);
      num3 = nurls.number(3);

      % Insert knots along the w direction
      if isempty(intins{3})
        coefs = nurls.coefs;
        knots{3} = nurls.knots{3};
        intervals{3}=nurls.intervals{3}; 
      else
        coefs = reshape(nurls.coefs,4*num1*num2,num3);
        [coefs,knots{3},intervals{3}] = nulintins(degree(3),coefs,nurls.knots{3},nurls.intervals{3},intins{3});
        num3 = size(coefs,2);
        coefs = reshape(coefs,[4 num1 num2 num3]);
      end

      % Insert knots along the v direction
      if isempty(intins{2})
        knots{2} = nurls.knots{2};
        intervals{2}=nurls.intervals{2}; 
      else
        coefs = permute(coefs,[1 2 4 3]);
        coefs = reshape(coefs,4*num1*num3,num2);
        [coefs,knots{2},intervals{2}] = nulintins(degree(2),coefs,nurls.knots{2},nurls.intervals{2},intins{2});
        num2 = size(coefs,2);
        coefs = reshape(coefs,[4 num1 num3 num2]);
        coefs = permute(coefs,[1 2 4 3]);
      end

      % Insert knots along the u direction
      if isempty(intins{1})
        knots{1} = nurls.knots{1};
        intervals{1}=nurls.intervals{1}; 
      else   
        coefs = permute(coefs,[1 3 4 2]);
        coefs = reshape(coefs,4*num2*num3,num1);
        [coefs,knots{1},intervals{1}] = nulintins(degree(1),coefs,nurls.knots{1},nurls.intervals{1},intins{1});
        coefs = reshape(coefs,[4 num2 num3 size(coefs,2)]);
        coefs = permute(coefs,[1 4 2 3]);
      end
     
    elseif size(nurls.knots,2)==2
      % NURL represents a surface
      num1 = nurls.number(1);
      num2 = nurls.number(2);

      % Insert knots along the v direction
      if isempty(intins{2})
        coefs = nurls.coefs;
        knots{2} = nurls.knots{2};
        intervals{2}=nurls.intervals{2}; 
      else
        coefs = reshape(nurls.coefs,4*num1,num2);
        [coefs,knots{2},intervals{2}] = nulintins(degree(2),coefs,nurls.knots{2},nurls.intervals{2},intins{2});
        num2 = size(coefs,2);
        coefs = reshape(coefs,[4 num1 num2]);
      end

      % Insert knots along the u direction
      if isempty(intins{1})
        knots{1} = nurls.knots{1};
        intervals{1}=nurls.intervals{1}; 
      else   
        coefs = permute(coefs,[1 3 2]);
        coefs = reshape(coefs,4*num2,num1);
        [coefs,knots{1},intervals{1}] = nulintins(degree(1),coefs,nurls.knots{1},nurls.intervals{1},intins{1});
        coefs = reshape(coefs,[4 num2 size(coefs,2)]);
        coefs = permute(coefs,[1 3 2]);
      end
     end
else
      % NURL represents a curve
      if isempty(intins)
        coefs = nurls.coefs;
        knots = nurls.knots;
        intervals=nurls.intervals;
      else
        [coefs, knots, intervals] = nulintins(degree, nurls.coefs, nurls.knots, nurls.intervals, intins);  
      end
end

% Construct new NURL
    coefs(1:3, :)=coefs(1:3, :)./repmat(coefs(4, :), 3, 1);
    inurl = nrlmake(coefs, knots, intervals, degree); 

end

%% demo
%  crv = nrbtestcrv;
%  crv=nrb2nrl(crv);
%  plot(crv.coefs(1,:),crv.coefs(2,:),'bo')
%  title('Interval insertion along test curve: curve and control points.');
%  hold on;
%  
%  nrlplot(crv,49);
%  
%  icrv = nrlintins(crv,[0.15 0.6 0.8] );
%  figure;
%  plot(icrv.coefs(1,:),icrv.coefs(2,:),'ro');
%  hold on;
%  nrlplot(crv,49);





