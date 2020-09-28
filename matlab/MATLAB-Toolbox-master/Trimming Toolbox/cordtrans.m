function [trans,npt] = cordtrans( pt,orig,theta )
% Get the local coordinates of given points in the local coordinate system
% from the global coordinate system. This function is only available for
% 2D, which means only rotation by z-axis. But the input pt should have 3D
% elements.

% Input:
%     pt: The 3D coordinates given points in the global cor sys, which is a
%       array or a matrix.(size(pt)=[numpts,3])
%     orig: The 3D coordinates of the local origin corresponding to the global cor sys.
%     theta: The rotation angle from the x-axis in global sys to the x-axis
%       in local sys, which is represented by rad.
% Output:
%     npt: The coordinates of the given point in the local cor sys.
%     trans: The trans matrix, which is a 4*4 matrix.

ttrans=vectrans(-orig); % from local sys to global sys
rtrans=vecrotz(-theta); % from local to global

pt=[pt';ones(1,length(pt))];
trans=rtrans*ttrans;
npt=trans*pt;
npt=npt';
npt(:,end)=[]; % Non-homogeneous coor.

end

