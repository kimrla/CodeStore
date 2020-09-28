function mt=cordtrans3D(orig,x,y,z)

% Get the local coordinates of given points in the local coordinate system
% from the global coordinate system in 3D space. This function is suitable
% for both righ-hand system and left-hand system.Notice that all the input
% are non-homogeneous coords, and the output is homogeneous matrix.

% Input:
%   orig: Origin point of new cord-sys.
%   x,y,z: Direction vector (Needn't to be unit) of 3 axes of the new cord-sys.
% Output:
%   mt: Transform matrix from global cord-sys to local sys, where
%       newpoints=mt*points, and size(points)=[4,numpoints].

x=x/norm(x);
y=y/norm(y);
z=z/norm(z);

T=[1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    -orig(:)',1];
T=T';
R=[x(:),y(:),z(:)];
R(:,end+1)=[0;0;0];
R(end+1,:)=[0,0,0,1];
R=R';

mt=R*T;





