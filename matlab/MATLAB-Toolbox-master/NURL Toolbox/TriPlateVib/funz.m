function [z, zx, zy, zxx, zxy, zyy]=funz(x, y)

% A test two dimensional function

z=sin(x).*sin(y);
zx=cos(x).*sin(y);
zy=sin(x).*cos(y);
zxx=-sin(x).*sin(y);
zyy=-sin(x).*sin(y);
zxy=cos(x).*cos(y);