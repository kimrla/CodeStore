function crv=nrlsrf2crv(srf, u, v)

% Extract a curve from a nurls surface
% 
% Calling Sequence:
%      crv=nrlsrfcrv(srf, u, [])
%      crv=nrlsrfcrv(srf, [], v)
%  
%  Input: 
%      u, v -  parametric coordinates
%
%  Output:
%      crv - a nurls curve
%

if isempty(u)
    u=srf.knots{1}; v=v*ones(size(u));
    knots=u; intervals=srf.intervals{1};
    order=srf.order(1);
elseif isempty(v)
    v=srf.knots{2}; u=u*ones(size(v));
    knots=v; intervals=srf.intervals{2};
    order=srf.order(2);
end
[ptw, w]=nrleval(srf, [u; v]);
pnts=ptw./repmat(w,3,1);
crv=nrlmake([pnts; w], knots, intervals, order);


%% Demo
% R=4;  u=[]; v=0.5;
% s1=0; s2=2*pi; t1=0; t2=pi; 
% center=[5, 5, 4]; 
% srf = nrlsphere(R, center, s1, s2, t1, t2); 
% figure; hold on; 
% nrlplot(srf, [100, 100]); 
% axis equal; view(3); 
% shading interp; 
% 
% % Extract a curve from a surface
% crv=nrlsrf2crv(srf, u, v);
% nrlplot(crv, 100);
% view(3); 



