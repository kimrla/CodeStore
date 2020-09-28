function varargout = nrlgdeval (nurl, tt)

% NRLGDEVAL: Evaluation of first and second derivatives of NURL by interpolation matrix. 
% 
% Calling Sequences:
%
%     [pnt, jac] = nrlgdeval (crv, tt)
%     [pnt, jac] = nrlgdeval (srf, {tu tv})
%     [pnt, jac] = nrlgdeval (vol, {tu tv tw})
%     [pnt, jac, hess] = nrlgdeval (crv, tt)
%     [pnt, jac, hess] = nrldgeval (srf, {tu tv})
%     [pnt, jac, hess] = nrldgeval (vol, {tu tv tw})
%
% INPUTS:
%
%   crv,   srf,   vol   - original NURL curve, surface or volume.
%   dcrv,  dsrf,  dvol  - NURL derivative representation of crv, srf 
%                          or vol (see nrbderiv2)
%   dcrv2, dsrf2, dvol2 - NURL second derivative representation of crv,
%                          srf or vol (see nrbderiv2)
%   tt     - parametric evaluation points
%            If the nurbs is a surface or a volume then tt is a cell
%            {tu, tv} or {tu, tv, tw} are the parametric coordinates
%
% OUTPUT:
%
%   pnt  - evaluated points.
%   jac  - evaluated first derivatives (Jacobian).
%   hess - evaluated second derivatives (Hessian).
%


if (nargin < 2)
  error ('nrlrdeval: wrong number of input parameters')
end

if (~isstruct(nurl))
  error('NURL representation is not structure!');
end

if (~strcmp(nurl.form,'L-NURL'))
  error('Not a recognised NURL representation');
end

% Check whether tt is a cell array with row vectors
tt=checktt(tt);

if (iscell(nurl.knots))
  if (size(nurl.knots,2) == 3)
  % NURL structure represents a volume
    
    if (nargout == 2)
        % first derivatives
        knots=nurl.knots; 
        order=nurl.order; 
        nt1 = length(tt{1}); 
        nt2 = length(tt{2}); 
        nt3 = length(tt{3}); 
        weights=squeeze(nurl.coefs(4, :, :, :)); 

        m=length(nurl.intervals{1}); 
        n=length(nurl.intervals{2}); 
        k=length(nurl.intervals{3}); 
        pnt=zeros(3, nt1, nt2, nt3); 
        dpu=pnt; dpv=pnt; dpw=pnt; 
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
                    [pnti, jaci] = nrlgintvdeval(order, weights(pp1, pp2, pp3), knotsijr, ttijr);
                    coefs=reshape(nurl.coefs(1:3, pp1, pp2, pp3), 3, []);

                    % The points
                    pti=reshape(coefs*pnti, 3, nt1i, nt2j, nt3r);
                    pnt(:, qq1, qq2, qq3)=pti(:, :, :, :); 

                    % First derivatives
                    dpui=reshape(coefs*jaci{1}, 3, nt1i, nt2j, nt3r);
                    dpu(:, qq1, qq2, qq3)=dpui(:, :, :, :); 
                    dpvi=reshape(coefs*jaci{2}, 3, nt1i, nt2j, nt3r);
                    dpv(:, qq1, qq2, qq3)=dpvi(:, :, :, :); 
                    dpwi=reshape(coefs*jaci{3}, 3, nt1i, nt2j, nt3r);
                    dpw(:, qq1, qq2, qq3)=dpwi(:, :, :, :); 
                end
            end
        end
        jac={dpu, dpv, dpw};
        
    elseif (nargout == 3)
        % second derivatives
        knots=nurl.knots; 
        order=nurl.order; 
        nt1 = length(tt{1}); 
        nt2 = length(tt{2}); 
        nt3 = length(tt{3}); 
        weights=squeeze(nurl.coefs(4, :, :, :)); 

        m=length(nurl.intervals{1}); 
        n=length(nurl.intervals{2}); 
        k=length(nurl.intervals{3}); 
        pnt=zeros(3, nt1, nt2, nt3); 
        dpu=pnt; dpv=pnt; dpw=pnt; 
        dpuu=pnt; dpuv=pnt; dpuw=pnt; 
        dpvv=pnt; dpvw=pnt; dpww=pnt; 
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
                    [pnti, jaci, hessi] = nrlgintvdeval(order, weights(pp1, pp2, pp3), knotsijr, ttijr);
                    coefs=reshape(nurl.coefs(1:3, pp1, pp2, pp3), 3, []);

                    % The points
                    pti=reshape(coefs*pnti, 3, nt1i, nt2j, nt3r);
                    pnt(:, qq1, qq2, qq3)=pti(:, :, :, :); 

                    % First derivatives
                    dpui=reshape(coefs*jaci{1}, 3, nt1i, nt2j, nt3r);
                    dpu(:, qq1, qq2, qq3)=dpui(:, :, :, :); 
                    dpvi=reshape(coefs*jaci{2}, 3, nt1i, nt2j, nt3r);
                    dpv(:, qq1, qq2, qq3)=dpvi(:, :, :, :); 
                    dpwi=reshape(coefs*jaci{3}, 3, nt1i, nt2j, nt3r);
                    dpw(:, qq1, qq2, qq3)=dpwi(:, :, :, :); 

                    % Second derivatives
                    dpuui=reshape(coefs*hessi{1,1}, 3, nt1i, nt2j, nt3r);
                    dpuu(:, qq1, qq2, qq3)=dpuui(:, :, :, :); 
                    dpuvi=reshape(coefs*hessi{1,2}, 3, nt1i, nt2j, nt3r);
                    dpuv(:, qq1, qq2, qq3)=dpuvi(:, :, :, :); 
                    dpuwi=reshape(coefs*hessi{1,3}, 3, nt1i, nt2j, nt3r);
                    dpuw(:, qq1, qq2, qq3)=dpuwi(:, :, :, :); 
                    dpvvi=reshape(coefs*hessi{2,2}, 3, nt1i, nt2j, nt3r);
                    dpvv(:, qq1, qq2, qq3)=dpvvi(:, :, :, :); 
                    dpvwi=reshape(coefs*hessi{2,3}, 3, nt1i, nt2j, nt3r);
                    dpvw(:, qq1, qq2, qq3)=dpvwi(:, :, :, :); 
                    dpwwi=reshape(coefs*hessi{3,3}, 3, nt1i, nt2j, nt3r);
                    dpww(:, qq1, qq2, qq3)=dpwwi(:, :, :, :); 
                end
            end
        end
        jac={dpu, dpv, dpw};
        hess={dpuu, dpuv, dpuw
                  dpuv, dpvv, dpvw
                  dpuw, dpvw, dpww};
    end

  elseif (size(nurl.knots,2) == 2)
    % NURL structure represents a surface 
    
    if (nargout == 2)
        % first derivatives
        knots=nurl.knots; 
        order=nurl.order; 
        nt1 = length(tt{1}); 
        nt2 = length(tt{2}); 
        weights=squeeze(nurl.coefs(4,:,:)); 

        m=length(nurl.intervals{1}); 
        n=length(nurl.intervals{2}); 
        pnt=zeros(3, nt1, nt2); 
        dpu=pnt; dpv=pnt;
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
                [Gi, jaci] = nrlgintvdeval(order, weights(pp1, pp2), knotsij, ttij);
                coefs=reshape(nurl.coefs(1:3, pp1, pp2), 3, []);
                
                % The points
                pti=reshape(coefs*Gi, 3, nt1i, nt2j);
                pnt(:, qq1, qq2)=pti(:, :, :); 
                
                % First derivatives
                dpui=reshape(coefs*jaci{1}, 3, nt1i, nt2j);
                dpu(:, qq1, qq2)=dpui(:, :, :); 
                dpvi=reshape(coefs*jaci{2}, 3, nt1i, nt2j);
                dpv(:, qq1, qq2)=dpvi(:, :, :); 
            end
        end
        jac={dpu, dpv};
        
    elseif (nargout == 3)
        % second derivatives
        knots=nurl.knots; 
        order=nurl.order; 
        nt1 = length(tt{1}); 
        nt2 = length(tt{2}); 
        weights=squeeze(nurl.coefs(4,:,:)); 

        m=length(nurl.intervals{1}); 
        n=length(nurl.intervals{2}); 
        pnt=zeros(3, nt1, nt2); 
        dpu=pnt; dpv=pnt;
        dpuu=pnt; dpvv=pnt; dpuv=pnt;
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
                [Gi, jaci, hessi] = nrlgintvdeval(order, weights(pp1, pp2), knotsij, ttij);
                coefs=reshape(nurl.coefs(1:3, pp1, pp2), 3, []);
                
                % The points
                pti=reshape(coefs*Gi, 3, nt1i, nt2j);
                pnt(:, qq1, qq2)=pti(:, :, :); 
                
                % First derivatives
                dpui=reshape(coefs*jaci{1}, 3, nt1i, nt2j);
                dpu(:, qq1, qq2)=dpui(:, :, :); 
                dpvi=reshape(coefs*jaci{2}, 3, nt1i, nt2j);
                dpv(:, qq1, qq2)=dpvi(:, :, :); 
                
                % Second derivatives
                dpuui=reshape(coefs*hessi{1,1}, 3, nt1i, nt2j);
                dpuu(:, qq1, qq2)=dpuui(:, :, :); 
                dpuvi=reshape(coefs*hessi{1,2}, 3, nt1i, nt2j);
                dpuv(:, qq1, qq2)=dpuvi(:, :, :); 
                dpvvi=reshape(coefs*hessi{2,2}, 3, nt1i, nt2j);
                dpvv(:, qq1, qq2)=dpvvi(:, :, :); 
            end
        end
        jac={dpu, dpv};
        hess={dpuu, dpuv
                   dpuv, dpvv};
    end
    
  end
else
  % NURL is a curve    
  
    if (nargout == 2)
        % first derivative
        m=length(nurl.intervals);
        pnts=cell(1,m-1); jac=pnts;
        for i=1:m-1
            a=nurl.intervals(i); b=nurl.intervals(i+1);
            pp=nurl.knots>=a &  nurl.knots<=b;
            qq=tt>=a &  tt<=b;
            [pti, jaci] = nrlgintvdeval(nurl.order, nurl.coefs(4,pp), nurl.knots(pp), tt(qq));
            pnts{i}=nurl.coefs(1:3,pp)*pti;
            jac{i}=nurl.coefs(1:3,pp)*jaci;
        end
        pnt=cell2mat(pnts);
        jac=cell2mat(jac);        
    elseif (nargout == 3)
        % second derivative
        m=length(nurl.intervals);
        pnts=cell(1,m-1); jac=pnts; hess=pnts;
        for i=1:m-1
            a=nurl.intervals(i); b=nurl.intervals(i+1);
            pp=nurl.knots>=a &  nurl.knots<=b;
            qq=tt>=a &  tt<=b;
            [pti, jaci, hessi] = nrlgintvdeval(nurl.order, nurl.coefs(4,pp), nurl.knots(pp), tt(qq));
            pnts{i}=nurl.coefs(1:3,pp)*pti;
            jac{i}=nurl.coefs(1:3,pp)*jaci;
            hess{i}=nurl.coefs(1:3,pp)*hessi;
        end
        pnt=cell2mat(pnts);
        jac=cell2mat(jac);        
        hess=cell2mat(hess);        
    end  
end

varargout{1} = pnt;
varargout{2} = jac;
if (nargout == 3)
  varargout{3} = hess;
end

end


%% Demo - curve
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
% 
% % Test nrlgdeval
% [p, dp, ddp]=nrlgdeval (crv, t);
% figure; 
% nrlplot(crv, 50); hold on;
% plot(p(1,:), p(2,:), 'ro');
% axis equal;
% 
% figure; hold on;
% t=linspace(0, 1, 51);
% [p1, dp1, ddp1]=nrlgdeval (crv, t);
% plot(dp1(1,:), dp1(2,:));
% quiver(dp(1,:), dp(2,:), ddp(1,:), ddp(2,:));
% axis equal; 

%% demo - surface
% R=1; N=6;
% s1=0; s2=pi; t1=0; t2=pi/2;
% center=[1, 1, 1];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% srf = nrlintins(srf, {0.5, 0.5});
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal; 
% hold on; 
% 
% % Test nrlgveval
% m=20; n=13; 
% s=linspace(0, 1, m); 
% t=linspace(0, 1, n); 
% p = nrlgeval(srf, {s, t}); 
% plot3(p(1,:), p(2,:), p(3,:), 'ro'); 
% 
% % Test nrlgdeval
% nurls=srf; tt={s, t};
% [pnt, jac, hess] = nrlgdeval (srf, tt);
% quiver3(pnt(1,:), pnt(2,:), pnt(3,:), jac{1}(1,:), jac{1}(2,:), jac{1}(3,:));
% quiver3(pnt(1,:), pnt(2,:), pnt(3,:), jac{2}(1,:), jac{2}(2,:), jac{2}(3,:));
% 
% m=40; n=30; 
% s=linspace(0, 1, m); 
% t=linspace(0, 1, n); 
% tt={s, t};
% [pnt1, jac1] = nrlgdeval (srf, tt);
% nt1=m; nt2=n;
% x=reshape(jac1{1}(1,:), nt1, nt2);
% y=reshape(jac1{1}(2,:), nt1, nt2);
% z=reshape(jac1{1}(3,:), nt1, nt2);
% figure;
% surf(x,y,z);
% hold on
% quiver3(jac{1}(1,:), jac{1}(2,:), jac{1}(3,:), hess{1,1}(1,:), hess{1,1}(2,:), hess{1,1}(3,:));
% quiver3(jac{1}(1,:), jac{1}(2,:), jac{1}(3,:), hess{1,2}(1,:), hess{1,2}(2,:), hess{1,2}(3,:));
% axis equal;

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
% 
% % Test nrlgdeval
% [pnt, jac] = nrlgdeval (vol, tt);
% figure; 
% nrlplot(vol, [20, 21, 8]);
% hold on;
% quiver3(pnt(1,:), pnt(2,:), pnt(3,:), jac{1}(1,:), jac{1}(2,:), jac{1}(3,:));
% quiver3(pnt(1,:), pnt(2,:), pnt(3,:), jac{2}(1,:), jac{2}(2,:), jac{2}(3,:));
% quiver3(pnt(1,:), pnt(2,:), pnt(3,:), jac{3}(1,:), jac{3}(2,:), jac{3}(3,:));
% axis equal;






