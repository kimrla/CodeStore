function [crv1, crv2]=bertrand(crv, dt, n)

% Get double curves of a curve with a distance (known as Bertrand curves)
% 
%  Inputs: 
%
%      crv - the curve to be transformed
%      dt - the distance between the two new curves
%       n - the number of knots to be inserted in each interval 
%            to improve smoothness (default value is 0)
% 
% Output: 
%
%      crv1, crv2 - the two new curves get by double
%  
%  Examples:  
%     crv=nrlcirc(2); dt=0.5;
%     nrlplot(crv, 100);
%     [crv1, crv2]=bertrand(crv, dt);
%

isl=isnrlline(crv);

if isl
    pnts = nrleval(crv, [0 1])';
    jac=pnts(1,:)-pnts(2,:); 
    nrm(1)=-jac(2); 
    nrm(2)=jac(1);
    nrm(3)=0;
    nrm=vecnorm(nrm);
    dr1=dt*nrm(1:2);
    dr2=-dt*nrm(1:2);
    crv1 = plcrvtranslat(crv, dr1);
    crv2 = plcrvtranslat(crv, dr2);
else
    if nargin==2
        n=0;
    end

    % Knots insertion
    iknt=n*ones(1, length(crv.intervals)-1);
    crv = nrlkntins(crv, iknt);

    % Get points and normal vectors
    [pntw, w]=nrleval(crv, crv.knots);
    [~, dp]=nrldeval(crv, crv.knots);
    dn=[-dp(2,:); dp(1,:); 0*dp(3,:)];

    % Translate the curve to get new curves
    dt=dt/2;
    dn=vecnorm(dn);
    wr=repmat(w, 3, 1);
    coefs1=(pntw+dt*dn.*wr)./wr;
    coefs2=(pntw-dt*dn.*wr)./wr;

    % Make curves
    crv1=nrlmake([coefs1; w], crv.knots, crv.intervals, crv.order);
    crv2=nrlmake([coefs2; w], crv.knots, crv.intervals, crv.order);
end

%% Demo - circle
% crv=nrlcirc(2); dt=0.5;
% nrlplot(crv, 100);
% 
% [crv1, crv2]=bertrand(crv, dt);
% 
% hold on;
% nrlplot(crv1, 100);
% nrlplot(crv2, 100);

%% Demo - splines
% p=2; dt=1; 
% x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
% y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
% crv=nrlspline(p, x, y); 
% nrlctrlplot(crv);
% 
% [crv1, crv2]=bertrand(crv, dt, 5);
% 
% hold on;
% nrlplot(crv1, 100);
% nrlplot(crv2, 100);




