function inurl = nrlkntins(nurl, iknts)
% 
% NRLKNTINS: Insert a single or multiple knots into a NURL curve,
%            surface or volume.
% 
% Calling Sequence:
% 
%   icrv = nrlkntins(crv,iukntvs);
%   isrf = nrlkntins(srf,{iukntvs ivkntvs});
%   ivol = nrlkntins(vol,{iukntvs ivkntvs iwkntvs});
% 
% INPUT:
% 
%   crv		: NURL curve, see nrlmake.
% 
%   srf		: NURL surface, see nrlmake.
% 
%   srf		: NURL volume, see nrlmake.
% 
%   iukntvs	: Vector of numbers of knots to be inserted along U direction.
% 
%   ivkntvs	: Vector of numbers of knotls to be inserted along V direction.
% 
%   iwkntvs	: Vector of numbers of knots to be inserted along W direction.
% 
% OUTPUT:
% 
%   icrv	: new NURL structure for a curve with knots inserted.
% 
%   isrf	: new NURL structure for a surface with knots inserted.
% 
%   ivol	: new NURL structure for a volume with knots inserted.
% 
% Description:
% 
%   Inserts knots into the NURL data structure. Knots along 
%   the V direction can only inserted into NURL surfaces, not 
%   curve that are always defined along the U direction. 
%   This function use the NUL function nulkntins. 
% 
% See also:
% 
%   nulintins
%


if nargin < 2
  error('Input argument must include the NURL and knots to be inserted');
end

if ~isstruct(nurl)
  error('NURL representation is not structure!');
end

if ~strcmp(nurl.form,'L-NURL')
  error('Not a recognised NURL representation');
end

degree = nurl.order;
nurl.coefs(1:3, :)=nurl.coefs(1:3, :).*repmat(nurl.coefs(4, :), 3, 1);

if iscell(nurl.knots)
     if size(nurl.knots,2)==3
      % NURL represents a volume
      num1 = nurl.number(1);
      num2 = nurl.number(2);
      num3 = nurl.number(3);

      % Insert knots along the w direction
      if max(iknts{3})<=0
        coefs = nurl.coefs;
        knots{3} = nurl.knots{3};
      else
        coefs = reshape(nurl.coefs,4*num1*num2,num3);
        [coefs,knots{3}] = nulkntins(degree(3),coefs,nurl.knots{3},nurl.intervals{3},iknts{3});
        num3 = size(coefs,2);
        coefs = reshape(coefs,[4 num1 num2 num3]);
      end

      % Insert knots along the v direction
      if max(iknts{2})<=0
        knots{2} = nurl.knots{2};
      else
        coefs = permute(coefs,[1 2 4 3]);
        coefs = reshape(coefs,4*num1*num3,num2);
        [coefs,knots{2}] = nulkntins(degree(2),coefs,nurl.knots{2},nurl.intervals{2},iknts{2});
        num2 = size(coefs,2);
        coefs = reshape(coefs,[4 num1 num3 num2]);
        coefs = permute(coefs,[1 2 4 3]);
      end

      % Insert knots along the u direction
      if max(iknts{1})<=0
        knots{1} = nurl.knots{1};
      else   
        coefs = permute(coefs,[1 3 4 2]);
        coefs = reshape(coefs,4*num2*num3,num1);
        [coefs,knots{1}] = nulkntins(degree(1),coefs,nurl.knots{1},nurl.intervals{1},iknts{1});
        coefs = reshape(coefs,[4 num2 num3 size(coefs,2)]);
        coefs = permute(coefs,[1 4 2 3]);
      end
     
    elseif size(nurl.knots,2)==2
      % NURL represents a surface
      num1 = nurl.number(1);
      num2 = nurl.number(2);

      % Insert knots along the v direction
      if max(iknts{2})<=0
        coefs = nurl.coefs;
        knots{2} = nurl.knots{2};
      else
        coefs = reshape(nurl.coefs,4*num1,num2);
        [coefs,knots{2}] = nulkntins(degree(2),coefs,nurl.knots{2},nurl.intervals{2},iknts{2});
        num2 = size(coefs,2);
        coefs = reshape(coefs,[4 num1 num2]);
      end

      % Insert knots along the u direction
      if max(iknts{1})<=0
        knots{1} = nurl.knots{1};
      else   
        coefs = permute(coefs,[1 3 2]);
        coefs = reshape(coefs,4*num2,num1);
        [coefs,knots{1}] = nulkntins(degree(1),coefs,nurl.knots{1},nurl.intervals{1},iknts{1});
        coefs = reshape(coefs,[4 num2 size(coefs,2)]);
        coefs = permute(coefs,[1 3 2]);
      end
     end
else
      % NURL represents a curve
      if max(iknts)<=0
        coefs = nurl.coefs;
        knots = nurl.knots;
      else
        [coefs, knots] = nulkntins(degree, nurl.coefs, nurl.knots, nurl.intervals, iknts);  
      end
end

% Construct new NURL
    coefs(1:3, :)=coefs(1:3, :)./repmat(coefs(4, :), 3, 1);
    inurl = nrlmake(coefs, knots, nurl.intervals, degree); 

end

%% demo
%   crv = nrbtestcrv;
%   crv=nrb2nrl(crv);
%   plot(crv.coefs(1,:),crv.coefs(2,:),'bo')
%   title('Knot insertion along test curve: curve and control points.');
%   hold on;
%   
%   nrlplot(crv,49);
%   
%   iknt=zeros(1, length(crv.intervals)-1);
%   iknt(1)=1; iknt(end)=1; 
%   icrv = nrlkntins(crv,iknt);
%   figure;
%   plot(icrv.coefs(1,:),icrv.coefs(2,:),'ro');
%   hold on;
%   nrlplot(crv,49);





