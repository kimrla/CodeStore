function nrlaxisplot(nrl,mn,c)

% NRLAXISPLOT: Plot the axises of a nurl geometry.
% 
% INPUT: 
%     nrb: NURL curve, surface or volume
%     mn: the number of Number of evaluation points, for a surface 
%            or volume, a row vector with the number of points along 
%            each direction. for example [20, 30]. Seen also: NRLPLOT
%       c:  Color of the two axes, for example 'cgm'
%
% EXAMPLE:  nrlaxisplot(crv, 10, 'c')
%                   nrlaxisplot(crv, 'c')
%                   nrlaxisplot(crv, 10)
%                   nrlaxisplot(crv)
%                   nrlaxisplot(srf, [10, 20], 'cg')
%                   nrlaxisplot(srf, [10, 20])
%                   nrlaxisplot(srf, 'cg')
%                   nrlaxisplot(vol, [10, 20, 15], 'cgm')
%                   nrlaxisplot(vol, [10, 20, 15])
%                   nrlaxisplot(vol, 'cgm')
%                   nrlaxisplot(vol)
% 

ii=numel(nrl.knots);

if nargin==1
    if iscell(nrl.knots) && ii==2
        mn=[100, 100];
    elseif iscell(nrl.knots) && ii==3
        mn=[100, 100, 100];
    else
        mn=100;
    end 
    c='cgm';
end

if nargin==2
    tt=mn;
    if ischar(tt)
        c=tt;
        if iscell(nrl.knots) && ii==2
            mn=[30, 30];
        elseif iscell(nrl.knots) && ii==3
            mn=[30, 30, 30];
        else
            mn=30;
        end 
    else
        c='cgm';
    end
end

if iscell(nrl.knots) && ii==2
    nrlplot(nrl,mn);
    hold on
    nrlctrlplot(nrl);
    [xyz, r] = nrlaxis(nrl);
    quiver3(xyz(1), xyz(2), xyz(3), r(1,1), r(2,1), r(3,1), c(1), 'LineWidth', 2.0);
    quiver3(xyz(1), xyz(2), xyz(3), r(1,2), r(2,2), r(3,2), c(2), 'LineWidth', 2.0);
elseif iscell(nrl.knots) && ii==3
    nrlplot(nrl,mn);
    hold on
    nrlctrlplot(nrl);
    [xyz, r] = nrlaxis(nrl);
    quiver3(xyz(1), xyz(2), xyz(3), r(1,1), r(2,1), r(3,1), c(1), 'LineWidth', 2.0);
    quiver3(xyz(1), xyz(2), xyz(3), r(1,2), r(2,2), r(3,2), c(2), 'LineWidth', 2.0);
    quiver3(xyz(1), xyz(2), xyz(3), r(1,3), r(2,3), r(3,3), c(3), 'LineWidth', 2.0);
else
    nrlplot(nrl,mn);
    hold on
    nrlctrlplot(nrl);
    [xyz, r] = nrlaxis(nrl);
    quiver3(xyz(1), xyz(2), xyz(3), r(1,1), r(2,1), r(3,1), c(1), 'LineWidth', 2.0);
end 


%% demo
% srf1 = nrl4surf([0.0 0.0 0.3],[1.0 0.0 -0.3],[0.0 1.0 -0.3],[1.0 1.0 0.3]);
% srf1 = nrlkntins(srf1,{0.5 0.5});
% crv1 =  nrlcirc(0.7,[0.5 -0.5 1.0],deg2rad(40),deg2rad(140));
% crv2 = nrlcirc(0.7,[0.5 0.5 1.0],deg2rad(40),deg2rad(140));
% srf2 = nrlruled (crv1,crv2);
% figure; nrlaxisplot(crv1);
% figure; nrlaxisplot(srf1);




