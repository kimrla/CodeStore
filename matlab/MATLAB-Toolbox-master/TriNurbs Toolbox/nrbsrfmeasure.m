function [Lu, Lv]=nrbsrfmeasure(srf, s, t)

% nrbsrfmeasure: Measure the length of parameric curves on a nurbs surface.
% 
% Calling Sequences:
% 
%     [Lu, Lv]=nrbsrfmeasure(srf)
% 
%     [Lu, Lv]=nrbsrfmeasure(srf, s)
% 
%     [Lu, Lv]=nrbsrfmeasure(srf, [], t)
% 
% INPUTS:
%
%      srf - A nurbs surface
%
% OUTPUT:
% 
%     Lu, Lv - The  length of parameric curves on the nurbs surface
% 

srf=nrb2nrl(srf);
if nargin==1
    s=srf.knots{1};
    t=srf.knots{2};
elseif nargin==2
    t=srf.knots{2};
elseif nargin==3 && isempty(s)
    s=srf.knots{1};
end
m=length(s); n=length(t);
Lu=zeros(1, m); Lv=zeros(1, n);    
for i=1:m
    crv=nrlsrfextract(srf, s(i), []);
    Lu(i)=nrlmeasure(crv);
end
for j=1:n
    crv=nrlsrfextract(srf, [], t(j));
    Lv(j)=nrlmeasure(crv);
end

%% demo
% % Create a nurbs surface and given node numbers (m, n)
% m=5; n=6;
% srf = nrbtestsrf; 
% s=linspace(0, 1, m);
% t=linspace(0, 1, n);
% 
% % Plot the two surfaces
% figure; hold on; 
% nrbplot(srf, [50, 50]);
% axis equal; view(3); 
% shading interp; 
% 
% % Measure the length of parameric curves on a nurbs surface
% [Lu, Lv]=nrbsrfmeasure(srf, s, t);





