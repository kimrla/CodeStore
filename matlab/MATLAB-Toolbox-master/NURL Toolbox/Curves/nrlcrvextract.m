function [p1, dp]=nrlcrvextract(crv)

% Extract the two end points of a curve and get 
%          the derivatives of the curve at the two points 
%
%  Inputs: 
%      crv - a nurl curve
%
%  Output: 
%      startendps - the strat and end points
%      derive - first order derivatives at the two points
%  
%  Examples: crv = nrlellip(2, 1, [0, 0], 0, 2*pi/3);
%           [startendpnts, derive]=nrlcrvextract(crv);
%           

[p1, dp] = nrldeval(crv, [0, 1]);


%% Demo
% crv = nrlellip(2, 1, [0, 0], 0, 2*pi/3);
% nrlplot(crv); hold on;
% [p1, dp] = nrldeval(crv, [0, 1]);
% quiver3(p1(1,:), p1(2,:), p1(3,:), dp(1,:), dp(2,:), dp(3,:));
% axis equal;




