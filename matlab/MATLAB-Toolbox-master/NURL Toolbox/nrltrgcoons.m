function srf = nrltrgcoons(crv1, crv2, crv3)
% 
% NRLCOONS: Construction of a triangle Coons patch.
% 
% Calling Sequence:
% 
%   srf = nrlcoons(crv1, crv2, crv3)
% 
% INPUT:
% 
%   crv1, crv2, crv3	: NURL curves defining the bottom U direction boundary of
% 		the constructed NURL surface.
%
% OUTPUT:
% 
%   srf		: Coons NURL surface patch.
% 
% Description:
% 
%   Construction of a bilinearly blended Coons surface patch from four NURL
%   curves that define the boundary.
%

if nargin ~= 3
  error('Incorrect number of input arguments');
end

nrls={crv1, crv2, crv3};

% Take nrls{1} as reference curve 1
tnrls = cell(1,4);
tnrls{1}=nrls{1};    
[xyz1, r1] = nrlaxis(nrls{1});
nrls(1) = [];
toler=(1e-3)*(norm(r1));

% Find curve 3 that have the same origin as curve 1
for j=1:2
    rnrlj=nrlreverse(nrls{j});
    xyz = nrlaxis(nrls{j});
    rxyz = nrlaxis(rnrlj);
    if norm(xyz-xyz1)<toler
        tnrls{3}=nrls{j};
        nrls(j) = []; 
        break;
    elseif norm(rxyz-xyz1)<toler
        tnrls{3}=rnrlj;
        nrls(j) = []; 
        break;
    end
end

% Find curve 2 whose origin is the end of curve 1
xyz1 = nrlaxis(nrlreverse(tnrls{1}));
rnrlj=nrlreverse(nrls{1});
xyz = nrlaxis(nrls{1});
if norm(xyz-xyz1)<toler
    tnrls{2}=nrls{1};
else
    tnrls{2}=rnrlj;
end

% Get the fourth curve to form coons surface
tnrls{4} = nrlline(crv1.coefs(1:3, 1)', crv1.coefs(1:3, 1)');
u1=tnrls{1}; u2=tnrls{3};
v1=tnrls{4}; v2=tnrls{2}; 

if (max (abs (nrleval (u1, u1.knots(1)) - nrleval (v1, v1.knots(1)))) > 1e-10 || ...
    max (abs (nrleval (u1, u1.knots(end)) - nrleval (v2, v2.knots(1)))) > 1e-10 || ...
    max (abs (nrleval (u2, u2.knots(1)) - nrleval (v1, v1.knots(end)))) > 1e-10 || ...
    max (abs (nrleval (u2, u2.knots(end)) - nrleval (v2, v2.knots(end)))) > 1e-10)
  error ('The four curves do not define a closed boundary')
end

r1 = nrlruled(u1, u2);
r2 = nrltransp(nrlruled(v1, v2));
t  = nrl4surf(u1.coefs(1:3,1), u1.coefs(1:3,end), u2.coefs(1:3,1), u2.coefs(1:3,end));

% raise all surfaces to a common degree
du = max([r1.order(1), r2.order(1), t.order(1)]);
dv = max([r1.order(2), r2.order(2), t.order(2)]);
r1 = nrldegelev(r1, [du - r1.order(1), dv - r1.order(2)]);
r2 = nrldegelev(r2, [du - r2.order(1), dv - r2.order(2)]);
t  = nrldegelev(t,  [du - t.order(1),  dv - t.order(2)]);

% merge the knot vectors, to obtain a common knot vector

% U intervals
k1 = r1.intervals{1};
k2 = r2.intervals{1};
k3 = t.intervals{1};
ku = unique([k1 k2 k3]);

% V intervals
k1 = r1.intervals{2};
k2 = r2.intervals{2};
k3 = t.intervals{2};
kv = unique([k1 k2 k3]);

% Insert intervals
r1 = nrlintins(r1, {ku, kv});
r2 = nrlintins(r2, {ku, kv});
t  = nrlintins(t,  {ku, kv});

% U knots
nu1=countknts(r1.intervals{1}, r1.knots{1}); 
nu2=countknts(r2.intervals{1}, r2.knots{1});
nu3=countknts(t.intervals{1}, t.knots{1});
nu=max([nu1; nu2; nu3]);
inu1=nu-nu1; inu2=nu-nu2; inu3=nu-nu3;

% V knots
nv1=countknts(r1.intervals{2}, r1.knots{2}); 
nv2=countknts(r2.intervals{2}, r2.knots{2});
nv3=countknts(t.intervals{2}, t.knots{2});
nv=max([nv1; nv2; nv3]);
inv1=nv-nv1; inv2=nv-nv2; inv3=nv-nv3;

% Insert knots
r1 = nrlkntins(r1, {inu1, inv1});
r2 = nrlkntins(r2, {inu2, inv2});
t  = nrlkntins(t,  {inu3, inv3});

% combine coefficient to construct Coons surface
coefs(1,:,:) = r1.coefs(1,:,:) + r2.coefs(1,:,:) - t.coefs(1,:,:);
coefs(2,:,:) = r1.coefs(2,:,:) + r2.coefs(2,:,:) - t.coefs(2,:,:);
coefs(3,:,:) = r1.coefs(3,:,:) + r2.coefs(3,:,:) - t.coefs(3,:,:);
coefs(4,:,:) = r1.coefs(4,:,:) + r2.coefs(4,:,:) - t.coefs(4,:,:);
srf = nrlmake(coefs, r1.knots, r1.intervals, r1.order);

%% Demo
% crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
% a=sin(pi/4);
% crv2 = nrlline([0.0 0.0 0.0]',[a a 0.0]');
% crv3 = nrlcirc(1, [0, 0], 0, pi/4);
% crv4=nrlline(crv1.coefs(1:3, 1)', crv1.coefs(1:3, 1)');
% 
% srf = nrltrgcoons(crv1, crv2, crv3);
% 
% nrlplot(srf, [20, 20]);
% view(2); axis equal;






