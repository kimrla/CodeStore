clear; clc;
close all;
% The mesh seed length (h0)
h0=1.1;

% Create a nurbs sphere
center=[2,2,4];
circ=nrbcirc(4, center, 0, pi);
srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
srf2=nrbtestsrf;

% Transform a nurbs surface into triangular representation
tsrf1=nrb2tri(srf1, h0); 
tsrf2=nrb2tri(srf2, h0); 

% Get the edges of a tri-nurbs surface that intersected with another tri-nurbs surface
p2t1=tnrbpts2tri(tsrf1);
p2t2=tnrbpts2tri(tsrf2);

% Get the intersection points of two tri-nurbs surfaces and sort them
[sed1, stri2, spts1, spts2, spnts1, spnts2]=tnrbintersects(tsrf1, tsrf2, p2t1, p2t2);
nr=length(sed1);

% Create parametric curves and planes
for r=1:nr
    crvs1{r}=nrbspline(spts1{r}, 3);
    crvs2{r}=nrbspline(spts2{r}, 3);
    tcrv1{r}=nrb2tri(crvs1{r}, tsrf1.seeds(2)); 
    tcrv2{r}=nrb2tri(crvs2{r}, tsrf2.seeds(2)); 
end

% The parametric plane of surface 2
tpra2=tnrbparasrf(tsrf2);

% Plot the results
figure; hold on; 
tnrbplot(tsrf1); 
tnrbplot(tsrf2); 
for r=1:nr
    plot3(spnts1{r}(:,1), spnts1{r}(:,2), spnts1{r}(:,3), 'ro'); 
    plot3(spnts2{r}(:,1), spnts2{r}(:,2), spnts2{r}(:,3), 'r*'); 
end
axis equal; view(3); 
title('Geometric grid'); 

figure; hold on;
triplot(tsrf1.delaunay, tsrf1.nodes(:,1), tsrf1.nodes(:,2));  
for r=1:nr    
    tnrbplot(tcrv1{r}); 
    plot(spts1{r}(:,1), spts1{r}(:,2), 'k.', 'MarkerSize', 13); 
%     plot(spts1{r}(:,1), spts1{r}(:,2), 'r', 'LineWidth', 1); 
end
title('Parametric mesh of surface 1');  
axis equal; 

figure; hold on; 
triplot(tsrf2.delaunay, tsrf2.nodes(:,1), tsrf2.nodes(:,2)); 
for r=1:nr
    tnrbplot(tcrv2{r}); 
    plot(spts2{r}(:,1), spts2{r}(:,2), 'k.', 'MarkerSize', 13); 
%     plot(spts2{r}(:,1), spts2{r}(:,2), 'r', 'LineWidth', 1); 
end
title('Parametric mesh of surface 2'); 
axis equal;

figure; hold on; 
tnrbplot(tpra2); 
title('The parametric domain of surface 2'); 
axis equal;




