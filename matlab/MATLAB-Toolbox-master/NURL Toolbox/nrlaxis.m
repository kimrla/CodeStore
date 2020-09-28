function [xyz, r] = nrlaxis(nrl)

% NRLAXIS: Get the axises of a nurl geometry.
% 
% INPUT: 
%     nrb: NURL curve, surface or volume
% 
% OUTPUT:
% 
%     xyz: axis origin of nurbs geometry
%     r:  vectors of U, V, W axes of nurbs geometry
%
% SEE ALSO:  NRLAXISPLOT
% 

ii=numel(nrl.knots);

if iscell(nrl.knots) && ii==2
    pnts = nrleval(nrl, {[0, 1], [0, 1]});
    xyz(:,1)=pnts(:,1,1);
    r(:,1)=pnts(:,2,1) - pnts(:,1,1);
    r(:,2)=pnts(:,1,2) - pnts(:,1,1);
elseif iscell(nrl.knots) && ii==3
    pnts = nrleval(nrl, {[0, 1], [0, 1], [0, 1]});
    xyz(:,1)=pnts(:,1,1,1);
    r(:,1)=pnts(:,2,1,1)-pnts(:,1,1,1);
    r(:,2)=pnts(:,1,2,1)-pnts(:,1,1,1);
    r(:,3)=pnts(:,1,1,2)-pnts(:,1,1,1);
else
    pnts=nrleval(nrl, [0,1]);
    xyz(:,1)=pnts(:,1);
    r(:,1)=pnts(:,2)-pnts(:,1);
end

%% demo
% crv1 =  nrlcirc(0.7,[0.5 -0.5 1.0],deg2rad(40),deg2rad(140));
% nrlplot(crv1,10);
% hold on
% nrlctrlplot(crv1);
% [xyz, r] = nrlaxis(crv1);
% quiver3(xyz(1), xyz(2), xyz(3), r(1,1), r(2,1), r(3,1), 'c', 'LineWidth', 2.0);


