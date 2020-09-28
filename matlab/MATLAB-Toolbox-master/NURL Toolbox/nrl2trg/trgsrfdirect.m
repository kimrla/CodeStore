function [srf, du]=trgsrfdirect(srf, dr)

% Rearrange the direction (axis) of a triangular nurls patch
% 
% Calling Sequences:
% 
%    srf=trgsrfdirect(srf)
%    srf=trgsrfdirect(srf, dr)
% 
% INPUT:
% 
%   srf		: NURLS surface, see nrlmake.
%
%   dr		: normal direction of the surface (default is [0, 0, 1]).
% 
% OUTPUT:
%
%   srf		:  rearranged NURLS surface  
%
%   du    :   an index of the shape  of the patch or the direction of u axis
%               du = -1 means the patch is an ellipse 
%               du = 0 means the patch is a quadrangle 
%               du = 1 or 2 is the axis NOT pointed to the overlapped point 
%

if nargin==1
    dr=[0, 0, 1];
end

% Find the overlapped point 
pnts = nrleval(srf, {[0, 1], [0, 1]});
[xyz, r] = nrlaxis(srf);
pnts=[pnts(:,1,1), pnts(:,2,1), pnts(:,2,2), pnts(:,1,2)];
[~, left]=FindUnique(pnts');
d0=[];
for i=1:3
    if ~isempty(left{i})
        d0=[d0, left{i}];
    end
end

if isempty(d0)
    % The patch is a quadrangle 
    [~, r] = nrlaxis(srf);
    dn=dot(cross(r(:,1)', r(:,2)'), dr);
    if dn<0
        srf = nrlpermute(srf, [2, 1]);
    end
    du=0;
    return;
elseif length(d0)==2
    % The patch is an ellipse 
    du=-1;
    return;
else
    % The patch is a triangle
    ctr=pnts(:, d0);
    dm=DistanceMatrix(xyz', pnts');
    tol=max(max(dm))*1e-6;
    dx=norm(r(:,1)); 
    dy=norm(r(:,2)); 
    di=min([dx, dy]);
    if di<tol
        srf=nrlreverse(srf);
    end
    [~, r] = nrlaxis(srf);
    dn=dot(cross(r(:,1)', r(:,2)'), dr);
    if dn<0
        srf = nrlpermute(srf, [2, 1]);
    end
    p2 = nrleval(srf, [1; 0]);
    if norm(p2(:)-ctr)<tol
        srf=nrlreverse(srf, 2);
        srf = nrlpermute(srf, [2, 1]);
    end
    du=1;
end


%% ! Demo
% % The number integration nodes (m, n)
% m=15; n=16;
% 
% % Creat a triangular nurls patch
% crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
% a=sin(pi/4);
% crv2 = nrlline([0.0 0.0 0.0]',[a a 0.0]');
% crv3 = nrlcirc(1, [0, 0], 0, pi/4);
% srf = nrltrgcoons(crv1, crv2, crv3);
% 
% % Rearrange its direction and plot
% [srf, du]=trgsrfdirect(srf);
% 
% figure; hold on;
% nrlplot(srf, [20, 20]);
% nrlaxisplot(srf);
% view(2); axis equal;
% 
% % Quiver plot
% figure; hold on;
% nrlplot(srf, [20, 20], 'ctrl');
% nrlplot(srf, [20, 20], 'quiver');


