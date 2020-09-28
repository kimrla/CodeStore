function nurl=nrb2nrl(nurbs)

% Transform nurbs curve, surface or volume into nurl
%  
%   Input:
%      nurbs - a nurbs curve, surface or volume
%
%  Output:
%     nurl -  a nurl curve, surface or volume
% 

if (~isstruct(nurbs))
  error('NURBS representation is not structure!');
end

if (~strcmp(nurbs.form,'B-NURBS'))
  error('Not a recognised NURBS representation');
end

if (iscell(nurbs.knots))
    % NURBS structure represents a volume or surface
    n=size(nurbs.knots,2);
    knots=cell(1,n); intervals=knots;
    for i=1:n
        [knots{i}, intervals{i}]=findknts(nurbs.knots{i}, nurbs.order(i));
        knots{i}(1)=0; knots{i}(end)=1;
        intervals{i}(1)=0; intervals{i}(end)=1;
    end
else
    % NURBS structure represents a curve
    [knots, intervals]=findknts(nurbs.knots, nurbs.order);
    knots(1)=0; knots(end)=1;
    intervals(1)=0; intervals(end)=1;
end
[pnts, w]=nrbeval(nurbs, knots);     
pnts = pnts./repmat(w, [3, 1]);
coefs=[pnts; w]; 
nurl=nrlmake(coefs, knots, intervals, nurbs.order-1); 



% Find knots and intervals from NURBS
function [iknots, intervals]=findknts(knots, nrbd)

u=FindUnique(knots');
m=length(u); s=0;
order=nrbd-1;
intervals=knots(u);
iknots=zeros(1, (m-1)*order+1);
for i=1:m-1
    a=intervals(i); b=intervals(i+1);
    kk=order+1; 
    tt=linspace(a, b, kk); 
    q=s+1:s+kk; 
    iknots(q)=tt(:); 
    s=s+kk-1; 
end


%% Demo - curve
% crvs = nrbtestcrv;
% crv=crvs;
% figure; nrbctrlplot(crvs);
% crvl=nrb2nrl(crv);
% figure; hold on;
% nrlplot(crvl, 1000, 'ctrl'); 
% axis equal;

%% Demo - surface
%  srf = nrbtestsrf;
%  nrbctrlplot(srf)
%  title('Test surface')
%  hold off
%  
%  srf = nrltestsrf;
%  figure; nrlplot(srf,[20 30],'ctrl')
%  title('Test surface')
%  hold off
%  axis equal;

%% Demo - volume
%  crv1 = nrbcirc (1, [0 0], 0, pi/2);
%  crv2 = nrbcirc (2, [0 0], 0, pi/2);
%  srf  = nrbruled (crv1, crv2);
%  vol  = nrbextrude (srf, [0 0 1]);
%  nrbplot (vol, [30 10 10])
%  title ('Extrusion of the quarter of a ring')
%  
%  vol=nrb2nrl(vol);
%  figure; nrblplot (vol, [30 10 10])
%  title ('Extrusion of the quarter of a ring')

%% demo - torque surface
% crvs = nrbtestcrv;
% crv=crvs;
% 
% figure; nrbctrlplot(crvs);
% 
% % Transform nurbs curve into nurls curve
% crvl=nrb2nrl(crv);
% 
% % Evalueate the curves and their derivatives
% figure; hold on;
% nrlplot(crvl, 1000, 'ctrl'); 
% axis equal;
% 
% % Load surfaces
% load torque srfs
% figure; hold on;
% for i=1:numel(srfs)
%     srfs(i).coefs(3,:,:)=0;
%     srfl(i)=nrb2nrl(srfs(i));
%     nrlplot(srfl(i), [10, 10],'ctrl');
% end
% axis equal;
% 
% for i=1:numel(srfl)
%     nrledgeplot(srfl(i), 10)
% end
% 
% % Extrude the surfaces
% figure; hold on;
% for i=1:numel(srfs)
%     vols(i)=nrlextrude(srfl(i), [0 0 0.1]);
%     nrlplot(vols(i), [10, 10, 10], 'ctrl');
%     nrledgeplot(vols(i), 10)
% end
% view(3); axis equal;




