function [i, j]=cmpoints(crv1,crv2)

% Get index of the common points of two curves
% 
%  Inputs: 
% 
%     crv1 - curve 1,
%     crv2 - curve 2
%  
%  Output: 
% 
%      i - the strat (1) or end point (2) number of curve 1
%      j - the strat (1) or end point (2) number of curve 2
% 
%  Examples: 
% 
%     crv1=nrlcirc(0.5,[0,0,0],pi/2,3*pi/2);
%     crv2=nrlline([0,0.5], [0.5,0.5]);
%     [i, j]=cmpoints(crv1,crv2);
%

stedpnts1=nrlcrvextract(crv1);
stedpnts2=nrlcrvextract(crv2);
DM = DistanceMatrix(stedpnts1', stedpnts2');
dist = nrlmeasure (crv1);
index=abs(DM)<dist*1e-3;
ii=[1,1;2,2]; jj=[1,2;1,2];
i=ii(index); j=jj(index);



%% Demo
% crv1=nrlcirc(0.5,[0,0,0],pi/2,3*pi/2);
% crv2=nrlline([0,0.5], [0.5,0.5]);
% [i, j]=cmpoints(crv1,crv2);
% nrlplot(crv1);
% hold on;
% nrlplot(crv2);




