function fh = nrlsrfplot(srf, mn, cmp, clr)

% NRBSRFPLOT: Plot a nurl surface and return its figure handle.
% 
% INPUT:
%
%    srf :  a nurl surface.
%    mn :  the number of grid points in x- (and y)-direction
%
% OUTPUT:
%
%    cmp :  colormap of the nurl surface.
%      clr :   color of the boundar curves
% 


% Draw the surface
if nargin ==2
    cmp = 'summer';
    clr = 'k';
end
m=mn(1); n=mn(2);
[x, y, z] = nrlxyz(srf, [m, n]);
colormap(cmp);
fh(1) = surf(x,y,z,'FaceColor','interp', 'EdgeColor','none');
hold on;
fh(5)=0;

% Draw the boundar curves
crvs=nrlextract(srf);
for j=1:2
    [x, y, z] = nrlxyz(crvs(j), mn(1));
    fh(j+1)=plot3(x, y, z, clr);
end
for j=1:2
    [x, y, z] = nrlxyz(crvs(j+2), mn(2));
    fh(j+3)=plot3(x, y, z, clr);
end


%% Demo
% R=1; N=6;
% s1=0; s2=pi/2; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;
% figure;
% fh = nrlsrfplot(srf,[15,15]);
% axis equal;






