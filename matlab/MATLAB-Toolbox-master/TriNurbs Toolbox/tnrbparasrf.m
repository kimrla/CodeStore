function tpra=tnrbparasrf(tsrf)

% tnrbparasrf: Transform the parametric domain of a tri-nurbs surface into a tri-nurbs surface
% 
% Calling Sequences:
% 
%     tpra=tnrbparasrf(tsrf)
% 
% INPUTS:
%
%      tsrf - Triangular representation of the nurbs (tri-nurbs) surface.
%
% OUTPUT:
% 
%     tpra - Triangular representation of the parametric domain 
%               of the nurbs (tri-nurbs) surface.
%

n=size(tsrf.nodes, 1);
points=zeros(n,3);
points(:,1:2)=tsrf.nodes;
coefs = cat(3,[0 0; 0 1],[1 1; 0 1]);
knots = {[0 0 1 1]  [0 0 1 1]};
plane=nrbmak(coefs,knots);
tpra.form='Tri-NURBS';
tpra.dim=2;
tpra.nurbs=plane;
tpra.numbers=tsrf.numbers;
tpra.seeds(1,1:3)=min(tsrf.seeds(2:end));
tpra.nodes=tsrf.nodes;
tpra.points=points;
tpra.delaunay=tsrf.delaunay;

%% demo
% % The mesh seed length (h0)
% h0=1.1;
% 
% % Create a tri-nurbs test surface
% srf=nrbtestsrf;
% tsrf=nrb2tri(srf, h0); 
% tpra=tnrbparasrf(tsrf);
% 
% % Plot the results
% figure; tnrbplot(tsrf); axis equal;
% figure; tnrbplot(tpra); 
% axis equal; view(2);


