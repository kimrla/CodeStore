function [p, w] = nrleval(nurl, tt, der)
% 
% NRLEVAL: Evaluate a NURL at parametric points.
% 
% Calling Sequences:
% 
%   [p,w] = nrleval(crv,ut,der)
%   [p,w] = nrleval(srf,{ut,vt},der)
%   [p,w] = nrleval(vol,{ut,vt,wt},der)
%   [p,w] = nrleval(srf,pts,der)
% 
% INPUT:
% 
%   crv		: NURL curve, see nrlmake.
% 
%   srf		: NURL surface, see nrlmake.
%
%   vol		: NURL volume, see nrlmake.
% 
%   ut		: Parametric evaluation points along U direction.
%
%   vt		: Parametric evaluation points along V direction.
% 
%   wt		: Parametric evaluation points along W direction.
%
%   pts     : Array of scattered points in parametric domain
%
%    der - orders of derivatives to be evaluated for NUL.
%             Note inclusion of this argument means NUL instead
% 		      of coordinates are evaluated and returned.
% 
% OUTPUT:
%
%   p		: Evaluated points on the NURL curve, surface or volume as 
% 		Cartesian coordinates (x,y,z). If w is included on the lhs argument
% 		list the points are returned as homogeneous coordinates (wx,wy,wz).
% 
%   w		: Weights of the homogeneous coordinates of the evaluated
% 		points. Note inclusion of this argument changes the type 
% 		of coordinates returned in p (see above).
% 
% Description:
% 
%   Evaluation of NURL curves, surfaces or volume at parametric points along  
%   the U, V and W directions. Either homogeneous coordinates are returned
%   if the weights are requested in the lhs arguments, or as Cartesian coordinates.
%   This function utilises nuldeval.
%  
% See also:
%  
%     nuldeval, nrlplot
%


if (nargin < 2)
  error('Not enough input arguments');
end

if nargin==2
    der=zeros(size(nurl.order));
end

foption = 1;    % output format 3D cartesian coordinates
if (nargout == 2)
  foption = 0;  % output format 4D homogenous coordinates 
end

if (~isstruct(nurl))
  error('NURL representation is not structure!');
end

if (~strcmp(nurl.form,'L-NURL'))
  error('Not a recognised NURL representation');
end

%% Transform NURL coefs into homogeneous coordinates (wx,wy,wz)
nurl.coefs(1:3, :)=nurl.coefs(1:3, :).*repmat(nurl.coefs(4, :), 3, 1);

%% Check whether tt is a cell array with row vectors
tt=checktt(tt);

if (iscell(nurl.knots))
  if (size(nurl.knots,2) == 3)
    %% NURL structure represents a volume

    num1 = nurl.number(1); 
    num2 = nurl.number(2); 
    num3 = nurl.number(3); 
    degree = nurl.order;     

    if (iscell(tt))
      nt1 = numel (tt{1});
      nt2 = numel (tt{2});
      nt3 = numel (tt{3});

      %% evaluate along the w direction
      val = reshape (nurl.coefs, 4*num1*num2, num3);
      val = nuldeval (degree(3), val, nurl.knots{3}, tt{3}, nurl.intervals{3}, der(3));
      val = reshape (val, [4 num1 num2 nt3]);

      %% Evaluate along the v direction
      val = permute (val, [1 2 4 3]);
      val = reshape (val, 4*num1*nt3, num2);
      val = nuldeval (degree(2), val, nurl.knots{2}, tt{2}, nurl.intervals{2}, der(2));
      val = reshape (val, [4 num1 nt3 nt2]);
      val = permute (val, [1 2 4 3]);

      %% Evaluate along the u direction
      val = permute (val, [1 3 4 2]);
      val = reshape (val, 4*nt2*nt3, num1);
      val = nuldeval (degree(1), val, nurl.knots{1}, tt{1}, nurl.intervals{1}, der(1));
      val = reshape (val, [4 nt2 nt3 nt1]);
      val = permute (val, [1 4 2 3]);
      pnts = val;

      p = pnts(1:3,:,:,:);
      w = pnts(4,:,:,:);
      if nargin==3
          if max(der)~=0
              return;
          end
      end
      if (foption)
        p = p./repmat(w,[3 1 1 1]);
      end

    else

      %% Evaluate at scattered points
      %% tt(1,:) represents the u direction
      %% tt(2,:) represents the v direction
      %% tt(3,:) represents the w direction

      st = size(tt);
      if (st(1) ~= 3 && st(2) == 3 && numel(st) == 2)
        tt = tt';
        st = size (tt);
      end
      nt = prod(st(2:end));

      tt = reshape (tt, [3, nt]);

      %% evaluate along the w direction
      val = reshape(nurl.coefs,4*num1*num2,num3);
      val = nuldeval(degree(3),val,nurl.knots{3},tt(3,:),nurl.intervals{3}, der(3));
      val = reshape(val,[4 num1 num2 nt]);

      %% evaluate along the v direction
      val2 = zeros(4*num1,nt);
      for v = 1:nt
        coefs = reshape(val(:,:,:,v),4*num1,num2);
        val2(:,v) = nuldeval(degree(2),coefs,nurl.knots{2},tt(2,v),nurl.intervals{2}, der(2));
      end
      val2 = reshape(val2,[4 num1 nt]);

      %% evaluate along the u direction
      pnts = zeros(4,nt);
      for v = 1:nt
        coefs = reshape (val2(:,:,v), [4 num1]);
        pnts(:,v) = nuldeval(degree(1),coefs,nurl.knots{1},tt(1,v),nurl.intervals{1}, der(1));
      end

      w = pnts(4,:);
      p = pnts(1:3,:);
      if nargin==3
          if max(der)~=0
              return;
          end
      end
      if (foption)
        p = p./repmat(w,[3, 1]);
      end

      if (numel(st) ~= 2)
        w = reshape (w, st(2:end));
        p = reshape (p, [3, st(2:end)]);
      end
    end

  elseif (size(nurl.knots,2) == 2)
    %% NURL structure represents a surface
  
    num1 = nurl.number(1);
    num2 = nurl.number(2);
    degree = nurl.order;

    if (iscell(tt))
      %% Evaluate over a [u,v] grid
      %% tt{1} represents the u direction
      %% tt{2} represents the v direction

      nt1 = length(tt{1});
      nt2 = length(tt{2});
    
      %% Evaluate along the v direction
      val = reshape(nurl.coefs,4*num1,num2);
      val = nuldeval(degree(2),val,nurl.knots{2},tt{2},nurl.intervals{2}, der(2));
      val = reshape(val,[4 num1 nt2]);
    
      %% Evaluate along the u direction
      val = permute(val,[1 3 2]);
      val = reshape(val,4*nt2,num1);
      val = nuldeval(degree(1),val,nurl.knots{1},tt{1},nurl.intervals{1}, der(1));
      val = reshape(val,[4 nt2 nt1]);
      val = permute(val,[1 3 2]);

      w = val(4,:,:);
      p = val(1:3,:,:);
      if nargin==3
          if max(der)~=0
              return;
          end
      end
      if (foption)
          p = p./repmat(w,[3 1 1]);
      end

    else

      %% Evaluate at scattered points
      %% tt(1,:) represents the u direction
      %% tt(2,:) represents the v direction

      st = size(tt);
      if (st(1) ~= 2 && st(2) == 2 && numel(st) == 2)
        tt = tt';
        st = size (tt);
      end
      nt = prod(st(2:end));

      tt = reshape (tt, [2, nt]);

      val = reshape(nurl.coefs,4*num1,num2);
      val = nuldeval(degree(2),val,nurl.knots{2},tt(2,:),nurl.intervals{2}, der(2));
      val = reshape(val,[4 num1 nt]);

      %% evaluate along the u direction
      pnts = zeros(4,nt);
      for v = 1:nt
          coefs = reshape (val(:,:,v), [4 num1]);
          pnts(:,v) = nuldeval(degree(1),coefs,nurl.knots{1},tt(1,v),nurl.intervals{1}, der(1));
      end

      w = pnts(4,:);
      p = pnts(1:3,:);
      if nargin==3
          if max(der)~=0
              return;
          end
      end
      if (foption)
          p = p./repmat(w,[3, 1]);
      end

      if (numel(st) ~= 2)
          w = reshape (w, st(2:end));
          p = reshape (p, [3, st(2:end)]);
      end
        
    end

  end
else

  %% NURL structure represents a curve
  %%  tt represent a vector of parametric points in the u direction

  if (iscell (tt) && numel (tt) == 1)
    tt = cell2mat (tt);
  end
  
  st = size (tt);
  
  val = nuldeval(nurl.order, nurl.coefs, nurl.knots, tt(:)', nurl.intervals, der);

  w = val(4,:);
  p = val(1:3,:);
  if nargin==3
      if max(der)~=0
          return;
      end
  end
  if foption
    p = p./repmat(w,3,1);
  end

  if (st(1) ~= 1 || numel(st) ~= 2)
    w = reshape (w, st);
    p = reshape (p, [3, st]);
  end

end

end


