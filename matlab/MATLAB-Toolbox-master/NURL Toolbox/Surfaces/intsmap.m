function [map, pits]=intsmap(tt, ut, pcrv)

% Get the map of trimmed surface
%
% Calling Sequence:
% 
%   [map, pits]=intsmap(tt, ut, pcrv)
% 
% INPUT:
% 
%    tt		: parametric points of the surface
% 
%    ut    : parametric points of the trimmed curve on thesurface
%
% OUTPUT:
% 
%   map  :  a map of the trimmed rigion (1 - boundary, 0 - inside)
%   pits   :  parametric points of the mapped (trim curve)
% 

v=ut; nt=length(ut);
m=length(tt{1}); n=length(tt{2});
dist=1/nt;
map=true(m,n); pits=[];
for k=1:m
    ps1=nrleval(pcrv, ut);
    u=tt{1}(k); 
    cc=ones(size(v));
    ps2=[u*cc; v; 0*cc];
    crv2 = nrlline([u, v(1), 0.0]',[u, v(end), 0.0]');
    [uu, d] =acrvints(ps1, ps2, ut);
    [~, pp, d]=lincrvsints(pcrv, crv2, dist, uu, d, 5);
    vk=squeeze(pp(2,1,:));
    for i=1:length(d)
        p=tt{2}>vk(i);
        map(k,p)=~map(k,p);
    end
    pits=[pits, squeeze(pp(:,1,:))];
end


%% Demo
% % Load trimmed surfaces
% load trisrfs tsrf1 tsrf2;
% srf=tsrf1.surface;
% crv=tsrf1.curve;
% pcrv=tsrf1.paras;
% 
% % First derivatives
% m=25; n=26; nt=50;
% t1=linspace(0, 1, m);
% t2=linspace(0, 1, n);
% ut=linspace(0, 1, nt);
% tt={t1, t2};
% [pnts, jac]=nrldeval(srf, tt);
% 
% % Plot parametric curves
% [v1, u1]=meshgrid(tt{2}, tt{1});
% figure; hold on;
% plot(u1, v1);
% plot(u1', v1');
% nrlplot(pcrv, 100);
% 
% % Find regions
% [map, pits]=intsmap(tt, ut, pcrv);
% plot(pits(1,:), pits(2,:), 'ro');
% 
% % Trimmed surface plot
% x=squeeze(pnts(1,:,:)); 
% y=squeeze(pnts(2,:,:)); 
% z=squeeze(pnts(3,:,:)); 
% x(~map)=nan;
% figure; surf(x, y, z);
% shading interp;






