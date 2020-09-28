function [force,C]=force_repel_bound(nrb,uvi,uvj,delta,mode,varargin)

% Calculate the repeling force between a pair of points in the parameter
% domain(1D or 2D) of a nurbs model. Notice that the output doesn't include
% the signal but just the magnitude of the force.

% Input:
% nrb: A NURBS structure of a curve or a surface.
% uvi: Par-cord of the objective point. If it's a surface, uvi is a 1*2 vector.
% uvj: Par-cord of the remaining point. If it's a surface, uvj is a 1*2 vector.
% delta: Coordinate step to calculate the local curvature. If it's a surface, delta is a 1*2 vector.
% mode: Determine whether to calculate the repeling force of 2 internal
%   points(mode='repel') or the attraction force of boundary points(mode='attract').
% k: Coefficient in the denominator which influences the efficient distance
%   of the neighbor points. k>=1 and the default is k=2.
% q: Charge of the two points. The default is 1.
% c: Normalization constant coefficient. The default is 1.
% Output:
% force: Force of the objective ith point. If it's a surface,force is a 1*2
%   vector of the components of u and v direction.
% C: Coefficient relating to the local curvature.

k=2;q=0.016;c=1;
if nargin==8
    k=varargin{1};
    q=varargin{2};
    c=varargin{3};
elseif nargin==7
    k=varargin{1};
    q=varargin{2};
elseif nargin==6
    k=varargin{1};
elseif nargin<5
    error('There are not enough input');
end

dim=length(nrb.number);
% pointi=nrbeval(nrb,uvi);
% pointj=nrbeval(nrb,uvj);
% rij=pointj-pointi;
rij=uvj-uvi; % Notice that the method uses the par-cords instead of the phy-cords.
% The default direction of the force is i->j, which means the force is
% attraction force from i to j.
h=norm(rij);
% force=zeros(1,dim);
% kp=zeros(1,dim);

% If the points are outside of the boundary after dealing with delta, then
% it needs to be compulsorily moved to the boundary.
uvi1=uvi+delta;uvi2=uvi-delta;
ll=logical(uvi1>1);
uvi1(ll)=1;
ll=logical(uvi2<0);
uvi2(ll)=0;
% Direction vector from pointi to piontj
ij=@(nrb,up1,up2)(nrbeval(nrb,up1)-nrbeval(nrb,up2));

if strcmp(mode,'repel')
    if dim==1 % NURBS curve
        cij=norm(ij(nrb,uvi1,uvi2));    
    else % NURBS surface
        % Respectively calculate the component coefficient of u and v direction
        cij(1)=norm(ij(nrb,[uvi1(1),uvi(2)],[uvi2(1),uvi(2)]));
        cij(2)=norm(ij(nrb,[uvi(1),uvi1(2)],[uvi(1),uvi2(2)]));       
    end
elseif strcmp(mode,'attract')
    % Calculate the coefficient Cwall in the formula
    cij=norm(ij(nrb,uvi,uvj));% NURBS curve or surface
    cij=cij/1000;% Set the attraction force very big
end        

ll=logical(cij>eps);   % ll(2)=0:logical vector;ll=[1,0]:NOT a logical vector.
if sum(ll)==0 % There are degraded points.// Consider the sphere surface.
    % Another method is to avoid the singular points before input the coordinates
    force=uvj*0;
    return;
else
    delta=delta(ll);cij=cij(ll);    
    kp=2*delta./cij;
    C=c*sqrt(kp);
    force=q*q*C.*rij(ll)/(h^(k+1));
end

% Direction of repeling force is opposite to that of the attraction force.
if strcmp(mode,'repel')
    force=-force;
end
end

