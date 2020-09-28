function [p, w] = nrlgeval(nurl, tt)

% NRLGEVAL: Evaluate a NURL at parametric points by interpolation matrix.
% 
% Calling Sequences:
% 
%   [p, w] = nrlgeval(crv, ut)
%   [p, w] = nrlgeval(srf, {ut,vt})
%   [p, w] = nrlgeval(vol, {ut,vt,wt})
% 
% INPUT:
% 
%   crv		: NURL curve, see nrbmak.
% 
%   srf		: NURL surface, see nrbmak.
%
%   vol		: NURL volume, see nrbmak.
% 
%   ut		: Parametric evaluation points along U direction.
%
%   vt		: Parametric evaluation points along V direction.
% 
%   wt		: Parametric evaluation points along W direction.
%
%   pts     : Array of scattered points in parametric domain
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

%% Check whether tt is a cell array with row vectors
tt=checktt(tt);

if (iscell(nurl.knots))
  if (size(nurl.knots,2) == 3)
    %% NURL structure represents a volume
    knots=nurl.knots; 
    order=nurl.order; 
    nt1 = length(tt{1}); 
    nt2 = length(tt{2}); 
    nt3 = length(tt{3}); 
    weights=squeeze(nurl.coefs(4, :, :, :)); 

    m=length(nurl.intervals{1}); 
    n=length(nurl.intervals{2}); 
    k=length(nurl.intervals{3}); 
    p=zeros(3, nt1, nt2, nt3); 
    w=zeros(1, nt1, nt2, nt3); 
    for i=1:m-1
        a1=nurl.intervals{1}(i); b1=nurl.intervals{1}(i+1); 
        pp1=knots{1}>=a1 & knots{1}<=b1;
        qq1=tt{1}>=a1 & tt{1}<=b1;
        nt1i=length(find(qq1)); 
        for j=1:n-1
            a2=nurl.intervals{2}(j); b2=nurl.intervals{2}(j+1);
            pp2=knots{2}>=a2 & knots{2}<=b2; 
            qq2=tt{2}>=a2 & tt{2}<=b2;  
            nt2j=length(find(qq2)); 
            for r=1:k-1
                a3=nurl.intervals{3}(r); b3=nurl.intervals{3}(r+1);
                pp3=knots{3}>=a3 & knots{3}<=b3; 
                qq3=tt{3}>=a3 & tt{3}<=b3;  
                nt3r=length(find(qq3)); 
                knotsijr={knots{1}(pp1), knots{2}(pp2), knots{3}(pp3)}; 
                ttijr={tt{1}(qq1), tt{2}(qq2), tt{3}(qq3)}; 
                [G, wi] = nrlgintveval(order, weights(pp1, pp2, pp3), knotsijr, ttijr);
                coefs=reshape(nurl.coefs(1:3, pp1, pp2, pp3), 3, []);
                pti=reshape(coefs*G, 3, nt1i, nt2j, nt3r);
                p(:, qq1, qq2, qq3)=pti(:, :, :, :); 
                w(1, qq1, qq2, qq3)=wi(:, :, :); 
            end
        end
    end
    
    if foption
            p = p./repmat(w,3,1);
    end
        
  elseif (size(nurl.knots,2) == 2)
    %% NURL structure represents a surface
    if (iscell(tt))
        knots=nurl.knots; 
        order=nurl.order; 
        nt1 = length(tt{1}); 
        nt2 = length(tt{2}); 
        weights=squeeze(nurl.coefs(4,:,:)); 

        m=length(nurl.intervals{1}); 
        n=length(nurl.intervals{2}); 
        p=zeros(3, nt1, nt2); 
        w=zeros(1, nt1, nt2); 
        for i=1:m-1
            a1=nurl.intervals{1}(i); b1=nurl.intervals{1}(i+1); 
            pp1=knots{1}>=a1 & knots{1}<=b1;
            qq1=tt{1}>=a1 & tt{1}<=b1;
            nt1i=length(find(qq1)); 
            for j=1:n-1
                a2=nurl.intervals{2}(j); b2=nurl.intervals{2}(j+1);
                pp2=knots{2}>=a2 & knots{2}<=b2; 
                qq2=tt{2}>=a2 & tt{2}<=b2;  
                nt2j=length(find(qq2)); 
                knotsij={knots{1}(pp1), knots{2}(pp2)}; 
                ttij={tt{1}(qq1), tt{2}(qq2)}; 
                [G, wi] = nrlgintveval(order, weights(pp1, pp2), knotsij, ttij);
                coefs=reshape(nurl.coefs(1:3, pp1, pp2), 3, []);
                pti=reshape(coefs*G, 3, nt1i, nt2j);
                p(:, qq1, qq2)=pti(:, :, :); 
                w(1, qq1, qq2)=wi(:, :); 
            end
        end
        
        if foption
            p = p./repmat(w,3,1);
        end
        
    end
  end
else

  %% NURL structure represents a curve
  %%  tt represent a vector of parametric points in the u direction

    m=length(nurl.intervals);
    pnts=cell(1,m-1); w=pnts;
    for i=1:m-1
        a=nurl.intervals(i); b=nurl.intervals(i+1);
        pp=nurl.knots>=a &  nurl.knots<=b;
        qq=tt>=a &  tt<=b;
        [G, w{i}] = nrlgintveval(nurl.order, nurl.coefs(4,pp), nurl.knots(pp), tt(qq));
        pnts{i}=nurl.coefs(1:3,pp)*G;
    end
    pnts=cell2mat(pnts);
    p=pnts; w=cell2mat(w);  
    
    if foption
        p = p./repmat(w,3,1);
    end
    
end

end

%% demo - curve
% % Major and minor semi-axes (a, b)
% % Start and end angles (sang, eang)
% % Coordinates of the center (center)
% % The number of knots (N)
% a=2; b=1; N=6;
% sang=0; eang=2*pi;
% center=[0, 0];
% 
% % Get elliptic arcs
% crv = nrlellip(a, b, center, sang, eang);
% nrlplot(crv, 50); hold on;
% nrlplot(crv, 11, 'quiver');
% axis equal; 
% 
% % Test nrlgeval
% n=11;
% t=linspace(0, 1, n);
% crv=nrlkntins(crv, [6, 6]);
% p = nrlgeval(crv, t);
% plot(p(1,:), p(2,:), 'ro');
% axis equal;

%% demo - surface
% R=1; N=6;
% s1=0; s2=2*pi; t1=0; t2=pi;
% center=[1, 1, 1];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% srf = nrlintins(srf, {0.5, 0.5});
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal; 
% hold on; 
% 
% % Test nrlgeval
% m=10; n=13; 
% s=linspace(0, 1, m); 
% t=linspace(0, 1, n); 
% p = nrlgeval(srf, {s, t}); 
% plot3(p(1,:), p(2,:), p(3,:), 'ro'); 

%% demo - volume
% R=1; N=6;
% s1=0; s2=pi/2; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;
% 
% vol  = nrlextrude (srf, [0.2 0.2 0.2]);
% vol=nrlintins(vol, {0.5, 0.5, 0.5});
% figure;
% nrlplot(vol, [20, 21, 8], 'ctrl');
% view(3); axis equal;
% 
% % Test nrlgeval
% nurls=vol; 
% m=10; n=13; k=4; 
% t1=linspace(0, 1, m); 
% t2=linspace(0, 1, n); 
% t3=linspace(0, 1, k); 
% tt={t1, t2, t3};
% p = nrlgeval(nurls, tt);
% 
% hold on;
% plot3(p(1,:), p(2,:), p(3,:), 'ro');


