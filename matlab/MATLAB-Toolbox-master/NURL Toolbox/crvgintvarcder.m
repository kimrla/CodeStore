function varargout = crvgintvarcder (order, coefs, knots, tt)

% crvgintvarcder: Evaluation of the first and second derivatives of NURL curve in arc length coordinates.
% 
% Calling Sequences:
%
%     [pnt, jac] = crvgintvarcder (order, coefs, knots, tt)
%     [pnt, jac, hess] = crvgintvarcder (order, coefs, knots, tt)
%
% INPUTS:
%
%      order - order of the nurls basis 
% 
%      coefs -  Cartesian coordinates and weights
% 
%      knots - knot vectors of an interval
% 
%      tt     - parametric evaluation points. 
%
% OUTPUT:
%
%   pnt  - interpolation matrix used to evaluate points.
%   jac  - interpolation matrix used to evaluate first derivatives (Jacobian).
%   hess - interpolation matrix used to evaluate second derivatives (Hessian).
%
%   Notice: arc length coordinates are used.
% 

if (iscell(knots) && numel(knots)>1)
    error('NURL structure is not a curve.');
else
    % NURL is a curve  
    weights=coefs(4,:);
    n=length(tt);
    
    % Get the weights of derivatives with respect to parametric coordinate
    [pnt, jaci] = nrlgintvdeval(order, weights, knots, tt);
    dpt=coefs(1:3,:)*jaci;
    
    % Get derivatives of arc length with respect to parametric coordinate
    ds=zeros(1,n); 
    for i=1:n
        ds(i)=norm(dpt(:,i));
    end

    % Get the weights of first derivative with respect to arc length
    jac=jaci./repmat(ds, [length(knots), 1]);
    
    % second derivative
    if (nargout == 3)
        % Get the weights of derivatives with respect to parametric coordinate
        [pnt, jaci, hessi] = nrlgintvdeval(order, weights, knots, tt);
        dpt=coefs(1:3,:)*jaci;
        ddpt=coefs(1:3,:)*hessi;
        
        % Get derivatives of arc length with respect to parametric coordinate
        ds=zeros(1,n); 
        for i=1:n
            ds(i)=norm(dpt(:,i));
        end
        dds=sum(dpt.*ddpt)./ds;

        % Get the weights of derivatives with respect to arc length
        jac=jaci./repmat(ds, [length(knots), 1]);
        q=jac.*repmat(dds, [length(knots), 1]);
        hess=(hessi-q)./repmat(ds.^2, [length(knots), 1]);        
    end
end

varargout{1} = pnt;
varargout{2} = jac;
if (nargout == 3)
    varargout{3} = hess;
end

end


%% Demo
% % Major and minor semi-axes (a, b)
% % Start and end angles (sang, eang)
% % Coordinates of the center (center)
% % The number of knots (N)
% a=2; b=1; N=6;
% sang=0; eang=pi;
% center=[0, 0];
% 
% % Get an elliptic arcs
% crv = nrlellip(a, b, center, sang, eang);
% figure; hold on;
% nrlplot(crv, 20);
% axis equal;
% 
% % Get curvature by parametric coordinates
% n=11;
% t=linspace(0, 1, n);
% [pnt0, cvt0]=curvature(crv, t);
% 
% % Test crvgintvarcder by solving curvature of the curve by acr length coordinates
% crv=nrlkntins(crv, 10);
% [pnti, jac, hess] = crvgintvarcder (crv.order, crv.coefs, crv.knots, t);
% cvt=zeros(1,n);
% ddpt=crv.coefs(1:3,:)*hess;
% for i=1:n
%     cvt(i)=norm(ddpt(:,i));
% end
% 
% % First derivative of the curve with respect to arc length
% pnt=crv.coefs(1:3,:)*pnti;
% dpt=crv.coefs(1:3,:)*jac;
% plot(pnt(1,:), pnt(2,:), 'ro');
% quiver(pnt(1,:), pnt(2,:), dpt(1,:), dpt(2,:));
% axis equal;
% 
% % Second derivative of the curve with respect to arc length
% t=linspace(0, 1, 30);
% [~, jac] = crvgintvarcder (crv.order, crv.coefs, crv.knots, t);
% dpt1=crv.coefs(1:3,:)*jac;
% figure; hold on;
% plot(dpt1(1,:), dpt1(2,:));
% quiver(dpt(1,:), dpt(2,:), ddpt(1,:), ddpt(2,:));
% axis equal;
% 
% % Exact curvature of the acr
% figure; hold on;
% t=linspace(0, 1, 11);
% s=t*(eang-sang)+sang;
% pte=[a*cos(s); b*sin(s)];
% cvte=((a^2*b^2)./(a^2*sin(s).^2 - b^2*sin(s).^2 + b^2).^3).^(1/2);
% plot3(pte(1,:), pte(2,:), cvte);
% plot3(pnt(1,:), pnt(2,:), cvt, 'ro');
%







