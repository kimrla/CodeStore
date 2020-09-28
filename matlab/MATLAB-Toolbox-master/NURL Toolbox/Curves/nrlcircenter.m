function center=nrlcircenter(crv)

% Get the center of a nurl circle
% 
%  Inputs: 
%      crv - a nurl curve of circle
%
%  Output: 
%     center - the center of the circle
%
%  Examples: 
%    center=[6;5.4;0]; 
%	 crv=nrlcirc(3, center', pi, 2*pi);
%    center=nrlcircenter(crv);
%

% Get the radius vectors of the curve
n=5; tt=linspace(0,1,n);
[pnts, dps] = nrldeval(crv,tt);
dps=vecnorm(dps);
pnt1=pnts(:,1); dr1=dps(:,1);
for i=2:n
    dr2=dps(:,i);
    if norm(cross(dr1, dr2))>0.02
        pnt2=pnts(:,i);
        break;        
    end
end

% Create two radius lines and get the center
ABC1 = planeline('pointnormal', pnt1(1:2)', dr1(1:2)');
ABC2 = planeline('pointnormal', pnt2(1:2)', dr2(1:2)');
center=planeline('intersection', ABC1, ABC2 );




