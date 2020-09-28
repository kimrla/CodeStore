function dr=neartangent(crv, pnt)

% Get tangent vector of an end of a line that is nearest to a point
%   
%  Inputs: 
%       crv - a nurls curve
%       pnt - a point
% 
%  Output: 
%       dr - tangent vector nearest to pnt
% 
%  Examples: pnt=[1;0;0]; 
%                       crv=nrlcirc(3, [0, 0, 0], 0, 2*pi/3);
%                       dr=neartangent(crv, pnt);
%

[pnts, dps]=nrlcrvextract(crv);
dm=DistanceMatrix(pnt', pnts');
[~, I]=min(dm);
switch I
    case 1
        dr=-dps(:,1);
    case 2
        dr=dps(:,2);
end


%% Demo
% pnt=[1;0;0]; 
% crv=nrlcirc(3, [0, 0, 0], 0, 2*pi/3);
% dr=neartangent(crv, pnt);
% figure; hold on;
% nrlplot(crv);
% quiver3(pnt(1), pnt(2), pnt(3), dr(1), dr(2), dr(3));


