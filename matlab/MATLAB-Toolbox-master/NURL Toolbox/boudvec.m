function [bp, tbp, nbp]=boudvec(srf)

% boudvec: boundary points and tangent and normal vectors
%
% Calling Sequences:
%
%     [bp, tbp, nbp]=boudvec(srf)
%
% INPUTS:
%
%      srf - a nurl surface
%
% OUTPUT:
%
%   bp  - boundary points
%   
%   tbp  - boundary tangent vectors
%
%   nbp  - boundary normal vectors
% 
%  Discriptions:
%      
%      This rutine is used for vibration anslysis of thin plate
%

[bp{1}, jac] = nrldeval (srf, [srf.knots{1}; zeros(1, srf.number(1))]);
tbp{1}=vecnorm(squeeze(jac{1}));
[bp{2}, jac] = nrldeval (srf, [srf.knots{1}; ones(1, srf.number(1))]);
tbp{2}=vecnorm(squeeze(jac{1}));
[bp{3}, jac] = nrldeval (srf, [zeros(1, srf.number(2)); srf.knots{2}]);
tbp{3}=vecnorm(squeeze(jac{2}));
[bp{4}, jac] = nrldeval (srf, [ones(1, srf.number(2)); srf.knots{2}]);
tbp{4}=vecnorm(squeeze(jac{2}));
nbp=cell(1, 4);
for i=1:4
    nbp{i}=[-tbp{i}(2,:); tbp{i}(1,:); tbp{i}(3,:)];
end


%% demo
% % Order of NURLS basis (order) and number of nodes (nodes)
% order=3*[1, 1]; nodes=[10, 6]; 
% 
% % Create a sector
% crv1 = nrlcirc(0.5, [0, 0], 0, pi/4);
% crv2 = nrlcirc(1, [0, 0], 0, pi/4);
% srf=nrlruled(crv1, crv2);
% srf = nrldegelev(srf, order-srf.order);
% iknts=nodes-srf.number;
% srf=nrlkntins(srf, {iknts(1), iknts(2)});
% figure; nrlplot(srf, [100, 100], 'ctrl');
% view(2); axis equal;
% 
% % Boundary points and tangent and normal vectors
% [bp, tbp, nbp]=boudvec(srf);
% 
% % Plot boundary points and tangent and normal vectors
% figure; hold on;
% nrlplot(srf, [30, 30]);
% shading interp;
% for i=1:4
%     plot(bp{i}(1,:), bp{i}(2,:), 'ro');
%     quiver(bp{i}(1,:), bp{i}(2,:), tbp{i}(1,:), tbp{i}(2,:));
%     quiver(bp{i}(1,:), bp{i}(2,:), nbp{i}(1,:), nbp{i}(2,:));
% end










