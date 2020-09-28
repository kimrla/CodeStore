function rx=vecrot3(point,vector,angle)

% 3D rotation transform of an arbitrary axis in Euler Space, comparing with
% vecrot.m, where the axis is only defined in original point.

% Input:
%   point: One ending point of the rotation axis.
%   vector: Direction vector of the rotation axis, point2-point1, which
%       needn't to be unit vector.
%   angle: Rotation angle. Anticlockwise is positive. Notice that the
%       situation that if the object is rotated by y-axis, the signal is
%       different from that of x,z-axis.
% Output:
%   rx: Rotation transform matrix, which is homogeneous matrix.

% The application of rx is: newpoint=rx*point, where point is a
% column-vector. Absolute vector.

vector=vector/norm(vector);
ax=vector(1);ay=vector(2);az=vector(3);
A=[ax*ax,ax*ay,ax*az;
    ay*ax,ay*ay,ay*az;
    az*ax,az*ay,az*az];
A_=[0,-az,ay;
    az,0,-ax;
    -ay,ax,0];
M=A+cos(angle)*(eye(3)-A)+sin(angle)*A_;
M=[M;0,0,0];M(:,end+1)=[0,0,0,1];

point=-point;
T = [1 0 0 point(1); 0 1 0 point(2); 0 0 1 point(3); 0 0 0 1];

rx=inv(T)*M*T;

end
