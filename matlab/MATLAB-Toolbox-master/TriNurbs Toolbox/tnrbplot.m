function tnrbplot(tnrb)

% tnrbplot: Plot a triangular representation of nurbs surface.
% 
% Calling Sequences:
% 
%     tnrbplot(tnrb)
% 
% INPUTS:
% 
%     tnrb - Triangular representation of the nurbs surface.
%

if (~isstruct(tnrb))
  error('Tri-NURBS representation is not structure!');
end

if (~strcmp(tnrb.form,'Tri-NURBS'))
  error('Not a recognised Tri-NURBS representation');
end

if (tnrb.dim == 2)
    %% Tri-NURBS structure represents a surface
    trisurf(tnrb.delaunay, tnrb.points(:,1), tnrb.points(:,2), tnrb.points(:,3));
elseif (tnrb.dim == 1)
    %% Tri-NURBS structure represents a curve
    plot3(tnrb.points(:,1), tnrb.points(:,2), tnrb.points(:,3), 'LineWidth', 1.2);
end


%% demo - curve and plane
% % The mesh seed length (h0)
% h0=0.5;
% 
% % Create a plane square and a plane cuve
% lin1=nrbline([0,1], [9,1]);
% lin2=nrbline([0,6], [9,6]);
% srf=nrbruled(lin1, lin2);
% crv=nrbtestcrv;
% 
% % Transform a nurbs surface into triangular representation
% tsrf=nrb2tri(srf, h0);
% tcrv=nrb2tri(crv, h0);
% 
% % Plot results
% figure; hold on; 
% tnrbplot(tsrf);
% tnrbplot(tcrv);
% axis equal;


%% demo - survace
% % The mesh seed length (h0)
% h0=0.8;
% 
% % Create a nurbs sphere
% circ=nrbcirc(4, [5,5,4], 0, pi);
% srf1=nrbrevolve(circ, [5,5,4], [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% 
% % Transform a nurbs surface into triangular representation
% tnrb1=nrb2tri(srf1, h0);
% tnrb2=nrb2tri(srf2, h0);
% 
% % Plot the results
% figure; hold on;
% tnrbplot(tnrb1);
% tnrbplot(tnrb2);
% axis equal; view(3);
% title('Geometric grid');




