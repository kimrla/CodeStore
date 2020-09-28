function nrlctrlplot (nurl)

% NRLCTRLPLOT: Plot a NURL entity along with its control points.
% 
% Calling Sequence:
% 
%   nlctrlplot (nurbs)
% 
% INPUT:
% 
%   nurbs: NURL curve, surface or volume, see nrbmak.
% 
% Example:
%
%   Plot the test curve and test surface with their control polygon and
%    control net, respectively
%
%   nrlctrlplot(nrltestcrv)
%   nrlctrlplot(nrltestsrf)
%

if (nargin < 1)
  error ('nrbctrlplot: Need a NURL to plot!');
end

% Default values
light='on';
cmap='summer';

colormap (cmap);

hold_flag = ishold;

if (iscell (nurl.knots))
  if (size (nurl.knots,2) == 3)
    nsub = 100;
    nrblplot (nurl, [nsub nsub nsub], 'light', light, 'colormap', cmap);
    hold on
% Plot the control points
    coefs = nurl.coefs(1:3,:,:,:);
    plot3 (coefs(1,:), coefs(2,:), coefs(3,:), 'r.','MarkerSize',20);

    % Plot the control net
    for ii = 1:size (nurl.coefs, 2)
         for jj = 1:size (nurl.coefs, 3)
              coefs = reshape (nurl.coefs(1:3,ii,jj,:), 3, []);
              plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'k--')
         end
         for kk = 1:size (nurl.coefs, 4)
              coefs = reshape (nurl.coefs(1:3,ii,:,kk), 3, []);
              plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'k--')
         end
    end
    for jj = 1:size (nurl.coefs, 3)
         for kk = 1:size (nurl.coefs, 4)
              coefs = reshape (nurl.coefs(1:3,:,jj,kk), 3, []);
              plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'k--')
         end
    end 
  elseif (size (nurl.knots,2) == 2) % plot a NURL surface

    nsub = 100;
    nrblplot (nurl, [nsub nsub], 'light', light, 'colormap', cmap);
    hold on

% And plot the control net
    for ii = 1:size (nurl.coefs, 2)
      coefs = reshape (nurl.coefs(1:3,ii,:), 3, []);
      plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'k--')
      plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'r.','MarkerSize',20)
    end
    for jj = 1:size (nurl.coefs, 3)
      coefs = reshape (nurl.coefs(1:3,:,jj), 3, []);
      plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'k--')
      plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'r.','MarkerSize',20)
    end
  end
else % plot a NURL curve
  nsub = 1000;
  nrblplot (nurl, nsub);
  hold on

% And plot the control polygon
  coefs = nurl.coefs(1:3,:);
  plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'k--')
  plot3 (coefs(1,:), coefs(2,:), coefs(3,:),'r.','MarkerSize',20)
end

if (~hold_flag)
  hold off
end

end

%% demo - curve
%  crv = nrltestcrv;
%  nrlctrlplot(crv)
%  title('Test curve')
%  hold off

%% demo - surface
%  srf = nrltestsrf;
%  nrlctrlplot(srf)
%  title('Test surface')
%  hold off

%% demo - comparison of nurbs and nurl curves
% % crv = nrbtestcrv;
% % crv=nrlcirc(1);
% % crv=nrl2nrb(crv);
% crv=nrbcirc(1);
% 
% figure; nrbctrlplot(crv);
% hold on;
% knots=RemDuplicate(crv.knots')';
% pnts=nrbeval(crv, knots);
% plot(pnts(1,:), pnts(2,:), 'kx', 'MarkerSize', 8, 'LineWidth', 0.8);
% 
% % Transform nurbs curve into nurls curve
% crvl=nrb2nrl(crv);
% 
% % Evalueate the curves and their derivatives
% figure; hold on;
% nrlplot(crvl, 1000, 'ctrl'); 
% axis equal;
% pnts=nrleval(crvl, crvl.intervals);
% plot(pnts(1,:), pnts(2,:), 'kx', 'MarkerSize', 8, 'LineWidth', 0.8);




