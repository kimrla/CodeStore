function srf = nrltestsrf

% NRBTESTSRF: Constructs a simple test surface.

% define a grid of control points
% in this case a regular grid of u,v points
% pnts(3,u,v)
%

srf=nrb2nrl(nrbtestsrf);

end

%% demo
% srf = nrltestsrf;
% nrlplot(srf,[20 30])
% title('Test surface')
% hold off

