function G = nuhgintveval(order, knots, tt, der)

% nuhgintveval: Evaluation nuh weights at parametric points.
% 
% Calling Sequences:
% 
%   G = nuhgintveval(order, knots, ut, der)
%   G = nuhgintveval(order, knots, {ut, vt}, der)
%   G = nuhgintveval(order, knots, {ut, vt, wt}, der)
%   G = nuhgintveval(order, knots, pts, der)
% 
% INPUT:
% 
%      order - order of the nurls basis 
% 
%      knots - knot vectors of an interval
%
%      tt     - parametric evaluation points. If the nurls is 
%            a surface or a volume then tt is a cell
%            {tu, tv} or {tu, tv, tw} are the parametric coordinates
%
%      der - orders of derivatives to be evaluated for NUH.
% 
% OUTPUT:
%
%   G		: Evaluated weights in Cartesian coordinates (x,y,z). 
%       If w is included on the lhs argument list the points are 
%       returned as homogeneous coordinates (wx,wy,wz).
% 
% Description:
% 
%   Evaluation of interpolation matrices or their derivatives for 
%   NUH curves, surfaces or volume at parametric points along  
%   the U, V and W directions. 
%  
% See also:
%  
%     nuhintvmat
%

if nargin==3
    der=zeros(size(order));
end

% Check whether tt is a cell array with row vectors
tt=checktt(tt);

if (iscell(knots))
  if (size(knots,2) == 3)
    %% NUH structure represents a volume

    num1 = length(knots{1});
    num2 = length(knots{2});
    num3 = length(knots{3});
    degree = order;     

    if (iscell(tt))      
      %% Evaluate weighting matrix for interpolation
      Gu=nuhintvmat(tt{1}, knots{1}, degree(1), der(1));
      Gv=nuhintvmat(tt{2}, knots{2}, degree(2), der(2));
      Gw=nuhintvmat(tt{3}, knots{3}, degree(3), der(3));
      G=kron(Gw, kron(Gv, Gu)); 
    else
        
      %% Evaluate at scattered points
      st = size(tt);
      if (st(1) ~= 3 && st(2) == 3 && numel(st) == 2)
        tt = tt';
        st = size (tt);
      end
      nt = prod(st(2:end));
      tt = reshape (tt, [3, nt]);
      
      %% Evaluate weighting matrix for interpolation
      Gu=nuhintvmat(tt(1,:), knots{1}, degree(1), der(1));
      Gv=nuhintvmat(tt(2,:), knots{2}, degree(2), der(2));
      Gw=nuhintvmat(tt(3,:), knots{3}, degree(3), der(3));
      G=zeros(num1*num2*num3, nt);
      for k=1:num3
          for j=1:num2
              for i=1:num1
                  p=(k-1)*num2*num1+(j-1)*num1+i; 
                  G(p,:)=Gu(i,:).*Gv(j,:).*Gw(k,:); 
              end
          end
      end
    end

  elseif (size(knots,2) == 2)
    %% NUH structure represents a surface
  
    num1 = length(knots{1});
    num2 = length(knots{2});
    degree = order;

    if (iscell(tt))
      %% Evaluate over a [u,v] grid      
      %% Evaluate weighting matrix for interpolation
      Gu=nuhintvmat(tt{1}, knots{1}, degree(1), der(1));
      Gv=nuhintvmat(tt{2}, knots{2}, degree(2), der(2));
      G=kron(Gv, Gu); 
    else
        %% Evaluate at scattered points

      st = size(tt);
      if (st(1) ~= 2 && st(2) == 2 && numel(st) == 2)
        tt = tt';
        st = size (tt);
      end
      nt = prod(st(2:end));
      tt = reshape (tt, [2, nt]);
      
      Gu=nuhintvmat(tt(1, :), knots{1}, degree(1), der(1));
      Gv=nuhintvmat(tt(2, :), knots{2}, degree(2), der(2));
      G=zeros(num1*num2, nt);
      for j=1:num2
          for i=1:num1
              p = (j-1)*num1+i;
              G(p,:) = Gu(i,:).*Gv(j,:);
          end
      end
    end
  end
else

  %% NUH basis represents a curve
  %%  tt represent a vector of parametric points in the u direction
  
  if (iscell (tt) && numel (tt) == 1)
      tt = cell2mat (tt);
  end
  
  G=nuhintvmat(tt, knots, order, der);

end

end


%% test of curve
% % A sine function and its derivatives
% fun = @(t) sin(t);
% dfun =@(t) cos(t);
% ddfun =@(t) -sin(t);
% dddfun =@(t) -cos(t);
% 
% % Oder of basis (order) and derivatives (der)
% order=2; der=0;
% 
% % knots (N), number of sampling points (m)
% n=10; m=21;
% 
% % Get knots (u) and nodes (t)
% u=linspace(-1, 1, n);
% t=linspace(-1, 1, m);
% 
% % Exact values
% y0=fun(u); dy0=dfun(u); 
% y=fun(t); dy=dfun(t); ddy=ddfun(t); d3y=dddfun(t); 
% 
% % Get derivative by weighting matrix
% G=nuhgintveval(order, u, t, der);
% D1=nuhgintveval(order, u, t, der+1);
% D2=nuhgintveval(order, u, t, der+2);
% D3=nuhgintveval(order, u, t, der+3);
% 
% n=length(u); 
% p=1:2:2*n-1; q=2:2:2*n; 
% yr=zeros(2*n,1);
% yr(p)=y0; yr(q)=dy0; 
% yv=G'*yr;
% ydv=D1'*yr;
% yd2v=D2'*yr;
% yd3v=D3'*yr;
% 
% % Plot results
% figure; hold on;
% plot(t, y); plot(t, yv, 'ro'); 
% title('Function values');
% 
% figure; plot(t, y-yv'); 
% title('Error of function values');
% 
% figure; hold on;
% plot(t, dy); plot(t, ydv, 'ro'); 
% title('First derivatives');
% 
% figure; plot(t, dy-ydv'); 
% title('Error of first derivatives');
% 
% figure; hold on;
% plot(t, ddy); plot(t, yd2v, 'ro'); 
% title('Second derivatives');
% 
% figure; plot(t, ddy-yd2v'); 
% title('Error of second derivatives');
% 
% figure; hold on;
% plot(t, d3y); plot(t, yd3v, 'ro'); 
% title('Third derivatives');
% 
% figure; plot(t, d3y-yd3v'); 
% title('Error of third derivatives');





