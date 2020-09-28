function varargout = nuhgintvdeval (order, knots, tt)

% nuhgintvdeval: Evaluation matrices of first and second derivatives of NURL curve, surface or volume.
% 
% Calling Sequences:
%
%     pnti = nuhgintvdeval (order, knots, tt)
%     pnti = nuhgintvdeval (order, knots, {tu tv})
%     pnti = nuhgintvdeval (order, knots, {tu tv tw})
%     [pnti, jaci] = nuhgintvdeval (order, knots, tt)
%     [pnti, jaci] = nuhgintvdeval (order, knots, {tu tv})
%     [pnti, jaci] = nuhgintvdeval (order, knots, {tu tv tw})
%     [pnti, jaci, hessi] = nuhgintvdeval (order, knots, tt)
%     [pnti, jaci, hessi] = nuhgintvdeval (order, knots, {tu tv})
%     [pnti, jaci, hessi] = nuhgintvdeval (order, knots, {tu tv tw})
%
% INPUTS:
%
%      order - order of the nurls basis 
% 
%      knots - knot vectors of an interval
% 
%      tt     - parametric evaluation points. If the nurls is 
%            a surface or a volume then tt is a cell
%            {tu, tv} or {tu, tv, tw} are the parametric coordinates
%
% OUTPUT:
%
%   pnti  - interpolation matrix used to evaluate points.
%   jaci  - interpolation matrix used to evaluate first derivatives (Jacobian).
%   hessi - interpolation matrix used to evaluate second derivatives (Hessian).
%
% See also:
%  
%     nuhgintveval
%


if (iscell(knots))
    % NUH structure represents a volume
	if (size(knots,2) == 3)
        % interpolation
        pnti=nuhgintveval(order, knots, tt, [0,0,0]);

        % first and second derivative
        if (nargout > 1)
            jaci{1}=nuhgintveval(order, knots, tt, [1, 0, 0]);
            jaci{2}=nuhgintveval(order, knots, tt, [0, 1, 0]);
            jaci{3}=nuhgintveval(order, knots, tt, [0, 0, 1]);
        end
        if (nargout >2)  
            hessi{1,1}=nuhgintveval(order, knots, tt, [2, 0, 0]);
            hessi{1,2}=nuhgintveval(order, knots, tt, [1, 1, 0]);
            hessi{1,3}=nuhgintveval(order, knots, tt, [1, 0, 1]);
            hessi{2,2}=nuhgintveval(order, knots, tt, [0, 2, 0]);
            hessi{2,3}=nuhgintveval(order, knots, tt, [0, 1, 1]);
            hessi{3,3}=nuhgintveval(order, knots, tt, [0, 0, 2]);
            hessi{2,1}=hessi{1,2};
            hessi{3,1}=hessi{1,3};
            hessi{3,2}=hessi{2,3};
        end
	% NUH structure represents a surface
	elseif (size(knots,2) == 2)
        % interpolation
        pnti=nuhgintveval(order, knots, tt, [0, 0]);

        % first and second derivative
        if (nargout > 1)
            jaci{1}=nuhgintveval(order, knots, tt, [1, 0]);
            jaci{2}=nuhgintveval(order, knots, tt, [0, 1]);
        end
        if (nargout >2)  
            hessi{1,1}=nuhgintveval(order, knots, tt, [2, 0]);
            hessi{1,2}=nuhgintveval(order, knots, tt, [1, 1]);
            hessi{2,2}=nuhgintveval(order, knots, tt, [0, 2]);
            hessi{2,1}=hessi{1,2};
        end
	end
else
    % NUH is a curve  
    % interpolation
    pnti=nuhgintveval(order, knots, tt, 0);

    % first and second derivative
    if (nargout > 1)
        jaci=nuhgintveval(order, knots, tt, 1);
    end
    if (nargout >2)  
        hessi=nuhgintveval(order, knots, tt, 2);
    end  
end

varargout{1} = pnti;
if (nargout >1)
    varargout{2} = jaci;
end
if (nargout >2)
  varargout{3} = hessi;
end


%% demo - curve
% % A sine function and its derivatives
% fun = @(t) sin(t);
% dfun =@(t) cos(t);
% ddfun =@(t) -sin(t);
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
% y=fun(t); dy=dfun(t); ddy=ddfun(t); 
% 
% % Get derivative by weighting matrix
% [G, D1, D2]=nuhgintvdeval(order, u, t);
% 
% [pt, jc]=nuhindex(n);
% yr=zeros(2*n,1);
% yr(pt)=y0; yr(jc)=dy0; 
% yv=G'*yr;
% ydv=D1'*yr;
% yd2v=D2'*yr;
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


%% demo - surface
% fun = @ (x, y) sin(x).*sin(y);
% funx = @ (x, y) cos(x).*sin(y);
% funy = @ (x, y) sin(x).*cos(y);
% funxx = @ (x, y) -sin(x).*sin(y);
% funyy = @ (x, y) -sin(x).*sin(y);
% funxy = @ (x, y) cos(x).*cos(y);
% 
% % Get parametric points and function values
% order=[2, 2]; m=10; n=12;
% M=m+3; N=n+3;
% s=linspace(0, pi/2, m);
% t=linspace(0, pi/2, n);
% x=linspace(0, pi/2, M);
% y=linspace(0, pi/2, N);
% [T, S]=meshgrid(t, s);
% [Y, X]=meshgrid(y, x);
% F=fun(S, T); Fx=funx(S, T); 
% Fy=funy(S, T); Fxy=funxy(S, T); 
% 
% % Get weighting matrices
% [pnti, jaci, hessi] = nuhgintvdeval (order, {s, t}, {x, y});
% G=pnti; Gx=jaci{1}; Gy=jaci{2}; 
% Gxx=hessi{1,1}; Gyy=hessi{2,2}; 
% Gxy=hessi{1,2}; 
% 
% % Get grid indexes
% [pt, jc, hs] = nuhindex(m, n);
% 
% % Prepare the coefficients vector
% ht=zeros(1, 4*m*n);
% ht(pt)=F(:);
% ht(jc{1})=Fx(:);
% ht(jc{2})=Fy(:);
% ht(hs)=Fxy(:);
% 
% % Interpolation
% Fi=reshape(ht*G, M, N);
% figure; surf(X, Y, Fi); title('Fi');
% figure; surf(X, Y, Fi-fun(X, Y)); title('error of Fi');
% 
% Fxi=reshape(ht*Gx, M, N);
% figure; surf(X, Y, Fxi); title('Fxi');
% figure; surf(X, Y, Fxi-funx(X, Y)); title('error of Fxi');
% 
% Fyi=reshape(ht*Gy, M, N);
% figure; surf(X, Y, Fyi); title('Fyi');
% figure; surf(X, Y, Fyi-funy(X, Y)); title('error of Fyi');
% 
% Fxxi=reshape(ht*Gxx, M, N);
% figure; surf(X, Y, Fxxi); title('Fxxi');
% figure; surf(X, Y, Fxxi-funxx(X, Y)); title('error of Fxxi');
% 
% Fyyi=reshape(ht*Gyy, M, N);
% figure; surf(X, Y, Fyyi); title('Fyyi');
% figure; surf(X, Y, Fyyi-funyy(X, Y)); title('error of Fyyi');
% 
% Fxyi=reshape(ht*Gxy, M, N);
% figure; surf(X, Y, Fxyi); title('Fxyi');
% figure; surf(X, Y, Fxyi-funxy(X, Y)); title('error of Fxyi');

%% demo - volume 
% fun = @ (x, y, z) sin(x).*sin(y).*sin(z); 
% funx = @ (x, y, z) cos(x).*sin(y).*sin(z); 
% funy = @ (x, y, z) sin(x).*cos(y).*sin(z); 
% funz = @ (x, y, z) sin(x).*sin(y).*cos(z); 
% funxx = @ (x, y, z) -sin(x).*sin(y).*sin(z); 
% funyy = @ (x, y, z) -sin(x).*sin(y).*sin(z); 
% funzz = @ (x, y, z) -sin(x).*sin(y).*sin(z); 
% funxy = @ (x, y, z) cos(x).*cos(y).*sin(z); 
% funyz = @ (x, y, z) sin(x).*cos(y).*cos(z); 
% funxz = @ (x, y, z) cos(x).*sin(y).*cos(z); 
% funxyz = @ (x, y, z) cos(x).*cos(y).*cos(z); 
% 
% % Get parametric points and function values
% order=[2, 2, 2]; m=5; n=6; k=7;
% M=m+3; N=n+3; K=k+3;
% r=linspace(0, pi/2, m);
% s=linspace(0, pi/2, n);
% t=linspace(0, pi/2, k);
% x=linspace(0, pi/2, M);
% y=linspace(0, pi/2, N);
% z=linspace(0, pi/2, K);
% [S, R, T]=meshgrid(s, r, t);
% [Y, X, Z]=meshgrid(y, x, z);
% F=fun(R, S, T); Fx=funx(R, S, T); Fy=funy(R, S, T); Fz=funz(R, S, T); 
% Fxx=funxx(R, S, T); Fyy=funyy(R, S, T); Fzz=funzz(R, S, T); 
% Fxy=funxy(R, S, T); Fyz=funyz(R, S, T); Fxz=funxz(R, S, T); 
% Fxyz=funxyz(R, S, T); 
% 
% % Get weighting matrices
% [pnti, jaci, hessi] = nuhgintvdeval (order, {r, s, t}, {x, y, z});
% G=pnti; Gx=jaci{1}; Gy=jaci{2}; Gz=jaci{3}; 
% Gxx=hessi{1,1}; Gyy=hessi{2,2}; Gzz=hessi{3,3}; 
% Gxy=hessi{1,2}; Gxz=hessi{1,3}; Gyz=hessi{2,3}; 
% 
% % Get grid indexes
% [pt, jc, hs, td] = nuhindex (m, n, k);
% 
% % Prepare the coefficients vector
% ht=zeros(1, 8*m*n*k); 
% ht(pt)=F(:); 
% ht(jc{1})=Fx(:); 
% ht(jc{2})=Fy(:); 
% ht(jc{3})=Fz(:); 
% ht(hs{1})=Fxy(:); 
% ht(hs{2})=Fxz(:); 
% ht(hs{3})=Fyz(:); 
% ht(td)=Fxyz(:); 
% 
% % Interpolation
% Fi=reshape(ht*G, M, N, K);
% Fe=fun(X, Y, Z);
% figure; volcplot(X, Y, Z, Fe); title('Fe');
% figure; volcplot(X, Y, Z, Fi-Fe); title('error of Fi');
% 
% Fxi=reshape(ht*Gx, M, N, K);
% Fxe=funx(X, Y, Z);
% figure; volcplot(X, Y, Z, Fxe); title('Fxe');
% figure; volcplot(X, Y, Z, Fxi-Fxe); title('error of Fxi');
% 
% Fyi=reshape(ht*Gy, M, N, K);
% Fye=funy(X, Y, Z);
% figure; volcplot(X, Y, Z, Fye); title('Fye');
% figure; volcplot(X, Y, Z, Fyi-Fye); title('error of Fyi');
% 
% Fzi=reshape(ht*Gz, M, N, K);
% Fze=funz(X, Y, Z);
% figure; volcplot(X, Y, Z, Fze); title('Fze');
% figure; volcplot(X, Y, Z, Fzi-Fze); title('error of Fzi');
% 
% Fxyi=reshape(ht*Gxy, M, N, K);
% Fxye=funxy(X, Y, Z);
% figure; volcplot(X, Y, Z, Fxye); title('Fxye');
% figure; volcplot(X, Y, Z, Fxyi-Fxye); title('error of Fxyi');
% 
% Fxzi=reshape(ht*Gxz, M, N, K);
% Fxze=funxz(X, Y, Z);
% figure; volcplot(X, Y, Z, Fxze); title('Fxze');
% figure; volcplot(X, Y, Z, Fxzi-Fxze); title('error of Fxzi');
% 
% Fyzi=reshape(ht*Gyz, M, N, K);
% Fyze=funyz(X, Y, Z);
% figure; volcplot(X, Y, Z, Fyze); title('Fyze');
% figure; volcplot(X, Y, Z, Fyzi-Fyze); title('error of Fyzi');
% 
% Fxxi=reshape(ht*Gxx, M, N, K);
% Fxxe=funxx(X, Y, Z);
% figure; volcplot(X, Y, Z, Fxxe); title('Fxxe');
% figure; volcplot(X, Y, Z, Fxxi-Fxxe); title('error of Fxxi');
% 
% Fyyi=reshape(ht*Gyy, M, N, K);
% Fyye=funyy(X, Y, Z);
% figure; volcplot(X, Y, Z, Fyye); title('Fyye');
% figure; volcplot(X, Y, Z, Fyyi-Fyye); title('error of Fyyi');
% 
% Fzzi=reshape(ht*Gzz, M, N, K);
% Fzze=funzz(X, Y, Z);
% figure; volcplot(X, Y, Z, Fzze); title('Fzze');
% figure; volcplot(X, Y, Z, Fzzi-Fzze); title('error of Fzzi');



