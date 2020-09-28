function force=force_resist(v,varargin)

% Calculate the medium resistance of a point. The formula is
% F=-k*h^m*v/h, where h=|v| and k is the damping coefficient and m is the
% order of the velocity(m>=2).

% Input:
%   v: Velocity of the point at the current time, which is a 1*2 vector.
%   k: Damping coefficient and the default is 2.
%   mm: Order of the velocity to calculate the damping force. The default is 2.
% Output:
%   force: Resistance force with the +- signal.

k=2;mm=2;
if nargin==3
    k=varargin{1};
    mm=varargin{2};
elseif nargin==2
    k=varargin{1};
end

force=-k*(norm(v))^(mm-1)*v;% Direction of the resistance force is opposite to that of velocity.
end


