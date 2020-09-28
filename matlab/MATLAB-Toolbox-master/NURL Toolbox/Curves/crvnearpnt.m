function [d, um, pntm]=crvnearpnt(crv, pnt, um)

% Get the nearest distance of a point to a curve
%
%  Calling sequence:
%    
%      [d, um, pntm]=crvnearpnt(crv, pnt)
% 
%      [d, um, pntm]=crvnearpnt(crv, pnt, um)
%
%  Inputs: 
%
%     crv - the curve to be used
% 
%     pnt - the point to be used
%
%     um - approximated parametric point on the curve
%
%  Output:
%    
%     d - the nearest distance
%
%     um - parametric point on the curve
%
%     pntm - x, y, z coordinates on the curve
%
%  Examples:  
%                 p=2; pnt=[6;6;0]; 
%                 x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
%                 y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
%                 crv=nrlspline(p, x, y);
%                 [d, um, pntm]=crvnearpnt(crv, pnt);

if nargin==2
    % Get the length of the curve
    dist= nrlmeasure (crv);

    % Get the maximum curvature of the curve
    n=5;
    for i=1:3
        ut=linspace(0, 1, n);
        [pnts, cvt]=curvature(crv, ut);
        n=fix(dist*max(cvt));
        if n<5
            n=5;
        elseif n>100
            n=100;
        end
    end

    % Get the approximated nearest points according to distance
    dm = DistanceMatrix(pnts',pnt'); 
    [~, I]=min(dm); 
    um=ut(I); 
end

% Solve the nearest points using Newton-Raphson's method
isl=isnrlline(crv, 5);
if isl
    [pntm, jac] = nrldeval (crv, um); 
    dr=pntm-pnt;
    f=dot(jac, dr);
    df=sum(jac.^2);
    um=um-f/df;
    pntm = nrleval (crv, um); 
else
    for i=1:5
        [pntm, jac, hess] = nrldeval (crv, um); 
        dr=pntm-pnt;
        f=dot(jac, dr);
        df=dot(hess, dr) + sum(jac.^2);
        um=um-f/df;
    end
end
d=norm(pntm-pnt);


%% Demo
% pnt=[6;6;0]; p=2;
% x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
% y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
% crv=nrlspline(p, x, y);
% nrlctrlplot(crv);
% hold on;
% plot(pnt(1), pnt(2), 'bo');
% [d, um, pntm]=crvnearpnt(crv, pnt);
% plot(pntm(1), pntm(2), 'ko');
% [pm, dpm]=nrldeval(crv, um);
% dpm=4*vecnorm(dpm);
% quiver(pm(1), pm(2), -dpm(2), dpm(1));


