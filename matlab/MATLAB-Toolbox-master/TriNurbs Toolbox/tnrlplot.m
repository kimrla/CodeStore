function tnrlplot(tnrl)

% tnrbplot: Plot a triangular representation of nurl surface.
% 
% Calling Sequences:
% 
%     tnrlplot(tnrl)
% 
% INPUTS:
% 
%     tnrl - Triangular representation of the nurl surface.
%

if (~isstruct(tnrl))
  error('Tri-NURL representation is not structure!');
end

if (~strcmp(tnrl.form,'Tri-NURL'))
  error('Not a recognised Tri-NURL representation');
end

if (tnrl.dim == 2)
    %% Tri-NURL structure represents a surface
    trisurf(tnrl.delaunay, tnrl.points(:,1), tnrl.points(:,2), tnrl.points(:,3));
elseif (tnrl.dim == 1)
    %% Tri-NURL structure represents a curve
    plot3(tnrl.points(:,1), tnrl.points(:,2), tnrl.points(:,3), 'LineWidth', 1.2);
end


%% demo - curve and plane
% % The mesh seed length (h0)
% h0=0.5;
% 
% % Create a plane square and a plane cuve
% lin1=nrlline([0,1], [9,1]);
% lin2=nrlline([0,6], [9,6]);
% srf=nrlruled(lin1, lin2);
% crv=nrltestcrv;
% 
% % Transform a nurbs surface into triangular representation
% tsrf=nrl2tri(srf, h0);
% tcrv=nrl2tri(crv, h0);
% 
% % Plot results
% figure; hold on; 
% tnrlplot(tsrf);
% tnrlplot(tcrv);
% axis equal;


%% demo - survace
% % The mesh seed length (h0)
% h0=0.8;
% 
% % Create a nurbs sphere
% circ=nrlcirc(4, [5,5,4], 0, pi);
% srf1=nrlrevolve(circ, [5,5,4], [1,0,0], 2*pi);
% srf2=nrltestsrf;
% 
% % Transform a nurbs surface into triangular representation
% tnrl1=nrl2tri(srf1, h0);
% tnrl2=nrl2tri(srf2, h0);
% 
% % Plot the results
% figure; hold on;
% tnrlplot(tnrl1);
% tnrlplot(tnrl2);
% axis equal; view(3);
% title('Geometric grid');




