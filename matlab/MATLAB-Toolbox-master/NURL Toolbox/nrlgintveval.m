function [G, w] = nrlgintveval(order, weights, knots, tt, der)

% nrlgintveval: Evaluation nurl weights at parametric points.
% 
% Calling Sequences:
% 
%   [G, w] = nrlgintveval(order, weights, knots, ut, der)
%   [G, w] = nrlgintveval(order, weights, knots, {ut, vt}, der)
%   [G, w] = nrlgintveval(order, weights, knots, {ut, vt, wt}, der)
%   [G, w] = nrlgintveval(order, weights, knots, pts, der)
% 
% INPUT:
% 
%      order - order of the nurls basis 
% 
%      weights -  Cartesian weights of an interval
% 
%      knots - knot vectors of an interval
%
%      tt     - parametric evaluation points. If the nurls is 
%            a surface or a volume then tt is a cell
%            {tu, tv} or {tu, tv, tw} are the parametric coordinates
%
%      der - orders of derivatives to be evaluated for NUL.
%             Note inclusion of this argument means NUL instead
% 		      of coordinates are evaluated and returned.
% 
% OUTPUT:
%
%   G		: Evaluated weights in Cartesian coordinates (x,y,z). 
%       If w is included on the lhs argument list the points are 
%       returned as homogeneous coordinates (wx,wy,wz).
% 
%   w		: Weights of the homogeneous coordinates.
%       Note inclusion of this argument changes the type 
% 		of coordinates returned in p (see above).
% 
% Description:
% 
%   Evaluation of interpolation matrices or their derivatives for 
%   NURL curves, surfaces or volume at parametric points along  
%   the U, V and W directions. Either homogeneous coordinates are returned
%   if the weights are requested in the lhs arguments, or as Cartesian coordinates.
%   This function utilises nuldeval and nulintvmat.
%  
% See also:
%  
%     nrlgeval
%

if nargin==4
    der=zeros(size(order));
end

foption = 1;    % output format cartesian coordinates
if (nargout == 2)
	foption = 0;  % output format homogenous coordinates 
end

% Check whether tt is a cell array with row vectors
tt=checktt(tt);

if (iscell(knots))
  if (size(knots,2) == 3)
    %% NURL structure represents a volume

    num1 = length(knots{1});
    num2 = length(knots{2});
    num3 = length(knots{3});
    degree = order;     

    if (iscell(tt))
      nt1 = numel (tt{1});
      nt2 = numel (tt{2});
      nt3 = numel (tt{3});
      intervals{1}=[min(knots{1}), max(knots{1})]; 
      intervals{2}=[min(knots{2}), max(knots{2})]; 
      intervals{3}=[min(knots{3}), max(knots{3})]; 

      %% evaluate along the w direction
      val = reshape (weights, num1*num2, num3);
      val = nuldeval (degree(3), val, knots{3}, tt{3}, intervals{3}, der(3));
      val = reshape (val, [num1 num2 nt3]);

      %% Evaluate along the v direction
      val = permute (val, [1 3 2]);
      val = reshape (val, num1*nt3, num2);
      val = nuldeval (degree(2), val, knots{2}, tt{2}, intervals{2}, der(2));
      val = reshape (val, [num1 nt3 nt2]);
      val = permute (val, [1 3 2]);

      %% Evaluate along the u direction
      val = permute (val, [2 3 1]);
      val = reshape (val, nt2*nt3, num1);
      val = nuldeval (degree(1), val, knots{1}, tt{1}, intervals{1}, der(1));
      val = reshape (val, [nt2 nt3 nt1]);
      val = permute (val, [3 1 2]);
      
      %% Evaluate weighting matrix for interpolation
      Gu=nulintvmat(tt{1}, knots{1}, degree(1), der(1));
      Gv=nulintvmat(tt{2}, knots{2}, degree(2), der(2));
      Gw=nulintvmat(tt{3}, knots{3}, degree(3), der(3));
      G=kron(Gw, kron(Gv, Gu)); 
      w = val; 
      G = G.*repmat(reshape(weights, [], 1), 1, nt1*nt2*nt3);  
 
      if max(der)~=0
           return;
      end
      if (foption)
          G=G./repmat(reshape(w, 1, []), [num1*num2*num3, 1]); 
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
      intervals{1}=[min(knots{1}), max(knots{1})]; 
      intervals{2}=[min(knots{2}), max(knots{2})]; 
      intervals{3}=[min(knots{3}), max(knots{3})]; 

      %% evaluate along the w direction
      val = reshape(weights, num1*num2, num3);
      val = nuldeval(degree(3), val, knots{3}, tt(3,:), intervals{3}, der(3));
      val = reshape(val,[1 num1 num2 nt]);

      %% evaluate along the v direction
      val2 = zeros(num1, nt);
      for v = 1:nt
        coefs = reshape(val(:,:,:,v), num1, num2);
        val2(:,v) = nuldeval(degree(2), coefs, knots{2}, tt(2,v), intervals{2}, der(2));
      end
      val2 = reshape(val2,[1 num1 nt]);

      %% evaluate along the u direction
      pnts = zeros(1,nt);
      for v = 1:nt
        coefs = reshape (val2(:,:,v), [1 num1]);
        pnts(:,v) = nuldeval(degree(1), coefs, knots{1}, tt(1,v), intervals{1}, der(1));
      end
      w = pnts(1,:);
      
      %% Evaluate weighting matrix for interpolation
      Gu=nulintvmat(tt(1,:), knots{1}, degree(1), der(1));
      Gv=nulintvmat(tt(2,:), knots{2}, degree(2), der(2));
      Gw=nulintvmat(tt(3,:), knots{3}, degree(3), der(3));
      G=zeros(num1*num2*num3, nt);
      for k=1:num3
          for j=1:num2
              for i=1:num1
                  p=(k-1)*num2*num1+(j-1)*num1+i; 
                  G(p,:)=Gu(i,:).*Gv(j,:).*Gw(k,:); 
              end
          end
      end
      G = G.*repmat(reshape(weights, [], 1), 1, nt);  
 
      if max(der)~=0
           return;
      end
      if (foption)
          G=G./repmat(reshape(w, 1, []), [num1*num2*num3, 1]); 
      end
    end

  elseif (size(knots,2) == 2)
    %% NURL structure represents a surface
  
    num1 = length(knots{1});
    num2 = length(knots{2});
    degree = order;

    if (iscell(tt))
      %% Evaluate over a [u,v] grid
      %% tt{1} represents the u direction
      %% tt{2} represents the v direction

      nt1 = length(tt{1});
      nt2 = length(tt{2});
      intervals={[min(knots{1}), max(knots{1})], [min(knots{2}), max(knots{2})]};
      
      %% Evaluate along the v direction
      val = weights; 
      val = nuldeval(degree(2), val, knots{2}, tt{2}, intervals{2}, der(2));
    
      %% Evaluate along the u direction
      val = permute(val,[2 1]);
      val = nuldeval(degree(1), val, knots{1}, tt{1}, intervals{1}, der(1));
      val = permute(val,[2 1]);     
      
      %% Evaluate weighting matrix for interpolation
      Gu=nulintvmat(tt{1}, knots{1}, degree(1), der(1));
      Gv=nulintvmat(tt{2}, knots{2}, degree(2), der(2));
      G=kron(Gv, Gu); 
      w = val; 
      G = G.*repmat(reshape(weights, [], 1), 1, nt1*nt2);  
 
      if max(der)~=0
           return;
      end
      if (foption)
          G=G./repmat(reshape(w, 1, []), [num1*num2, 1]);
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

      %% Evaluate along the v direction
      intervals={[min(knots{1}), max(knots{1})], [min(knots{2}), max(knots{2})]}; 
      val = weights; 
      val = nuldeval(degree(2), val, knots{2}, tt(2, :), intervals{2}, der(2)); 
      val = reshape(val, [1 num1 nt]);

      %% evaluate along the u direction
      pnts = zeros(1, nt);
      for v = 1 : nt
          coefs = reshape (val(:, :, v), [1 num1]);
          pnts(:, v) = nuldeval(degree(1), coefs, knots{1}, tt(1,v), intervals{1}, der(1));
      end
      w = pnts;
      if (numel(st) ~= 2)
          w = reshape (w, st(2 : end));
      end
      
      Gu=nulintvmat(tt(1, :), knots{1}, degree(1), der(1));
      Gv=nulintvmat(tt(2, :), knots{2}, degree(2), der(2));
      G=zeros(num1*num2, nt);
      for j=1:num2
          for i=1:num1
              p = (j-1)*num1+i;
              G(p,:) = Gu(i,:).*Gv(j,:);
          end
      end
      G = G.*repmat(reshape(weights, [], 1), 1, nt);  
      
      if max(der)~=0
          return;
      end
      if (foption)
          G=G./repmat(reshape(w, 1, []), [num1*num2, 1]);
      end
    end
  end
else

  %% NURL basis represents a curve
  %%  tt represent a vector of parametric points in the u direction
  
  if (iscell (tt) && numel (tt) == 1)
      tt = cell2mat (tt);
  end
  
  intervals = [min(knots), max(knots)];
  w = nuldeval(order, weights, knots, tt, intervals, der);
  G=nulintvmat(tt, knots, order, der);
  
  G = G.*repmat(weights', 1, length(tt));  
  
  if max(der)~=0
      return;
  end
  if foption
      G=G./repmat(w,[length(knots), 1]);
  end

end

end


%% Demo - curve
% % Major and minor semi-axes (a, b)
% % Start and end angles (sang, eang)
% % Coordinates of the center (center)
% % The number of knots (N)
% a=2; b=1; N=6;
% sang=0; eang=pi;
% center=[0, 0];
% 
% % Get elliptic arcs
% crv = nrlellip(a, b, center, sang, eang);
% nrlplot(crv, 30); hold on;
% nrlplot(crv, 10, 'quiver');
% axis equal; 
% 
% % Test nrlgintveval
% n=10;
% t=linspace(0, 1, n);
% crv=nrlkntins(crv, 6);
% G = nrlgintveval(crv.order, crv.coefs(4,:), crv.knots, t);
% p=crv.coefs(1:3,:)*G;
% plot(p(1,:), p(2,:), 'ro');

%% demo - surface
% % Creat a sphere
% R=1; N=6;
% s1=0; s2=pi/2; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% srf=nrlkntins(srf, {2, 1});
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;
% 
% % Test nrlgintveval
% m=10; n=13;
% s=linspace(0, 1, m);
% t=linspace(0, 1, n);
% tt={s, t};
% knots=srf.knots;
% weights=squeeze(srf.coefs(4,:,:));
% der=[0, 0];
% order=srf.order;
% 
% num1 = length(knots{1});
% num2 = length(knots{2});
% degree = order;
% 
% nt1 = length(tt{1});
% nt2 = length(tt{2});
% 
% G = nrlgintveval(order, weights, knots, tt, der);
% coefs=reshape(srf.coefs(1:3,:), 3, []);
% p=coefs*G;
% hold on;
% plot3(p(1,:), p(2,:), p(3,:), 'ro');
% axis equal;
% 
% x=reshape(p(1,:), nt1, nt2);
% y=reshape(p(2,:), nt1, nt2);
% z=reshape(p(3,:), nt1, nt2);
% figure;
% nrlplot(srf, [30, 30], 'ctrl');
% hold on
% surf(x,y,z);
% axis equal;
% 
% [u, v]=meshgrid(s, t);
% pts=[u(:)'; v(:)'];
% G = nrlgintveval(order, weights, knots, pts, der);
% p=coefs*G;
% figure;
% nrlplot(srf, [100, 100], 'ctrl');
% hold on;
% plot3(p(1,:), p(2,:), p(3,:), 'ro');
% axis equal;

%% demo - volume
% % Creat a sphere
% R=1; N=6;
% s1=0; s2=pi/2; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% vol  = nrlextrude (srf, [0.2 0.2 0.2]);
% vol=nrlkntins(vol, {2, 1, 1});
% figure;
% nrlplot(vol, [20, 21, 8], 'ctrl');
% view(3); axis equal;
% 
% % Test nrlgintveval
% m=10; n=13; k=4;
% t1=linspace(0, 1, m);
% t2=linspace(0, 1, n);
% t3=linspace(0, 1, k);
% tt={t1, t2, t3};
% knots=vol.knots;
% weights=squeeze(vol.coefs(4,:,:,:));
% der=[0, 0, 0];
% order=vol.order;
% 
% num1 = length(knots{1});
% num2 = length(knots{2});
% num3 = length(knots{3});
% degree = order;
% 
% nt1 = length(tt{1});
% nt2 = length(tt{2});
% nt3 = length(tt{3});
% 
% % Get and plot results
% G = nrlgintveval(order, weights, knots, tt, der);
% coefs=reshape(vol.coefs(1:3,:), 3, []); 
% p=coefs*G; 
% figure; hold on; 
% nrblplot(vol, [20 20 8], 'light', 'on'); 
% plot3(p(1,:), p(2,:), p(3,:), 'ro'); 
% view(3); axis equal; 
% 
% X=reshape(p(1,:), nt1, nt2, nt3);
% Y=reshape(p(2,:), nt1, nt2, nt3);
% Z=reshape(p(3,:), nt1, nt2, nt3);
% figure;
% nrlplot(vol, [20, 20, 6], 'ctrl');
% hold on
% m=nt1; n=nt2; k=nt3;
% surf(squeeze(X(1,:,:)), squeeze(Y(1,:,:)), squeeze(Z(1,:,:))); 
% surf(squeeze(X(m,:,:)), squeeze(Y(m,:,:)), squeeze(Z(m,:,:))); 
% surf(squeeze(X(:,1,:)), squeeze(Y(:,1,:)), squeeze(Z(:,1,:))); 
% surf(squeeze(X(:,n,:)), squeeze(Y(:,n,:)), squeeze(Z(:,n,:))); 
% surf(squeeze(X(:,:,1)), squeeze(Y(:,:,1)), squeeze(Z(:,:,1))); 
% surf(squeeze(X(:,:,k)), squeeze(Y(:,:,k)), squeeze(Z(:,:,k))); 
% axis equal;
% 
% [u, v, w]=meshgrid(t1, t2, t3);
% pts=[u(:)'; v(:)'; w(:)'];
% G = nrlgintveval(order, weights, knots, pts, der);
% coefs=reshape(vol.coefs(1:3,:), 3, []); 
% p=coefs*G; 
% figure; hold on; 
% nrblplot(vol, [20 20 8], 'light', 'on'); 
% plot3(p(1,:), p(2,:), p(3,:), 'ro'); 
% view(3); axis equal; 








