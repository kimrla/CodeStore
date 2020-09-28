function plotedge(x, y, z, c)

% Plot the edge of a surface
% 
% Calling Sequence:
%
%      plotedge(x, y, z)
%      plotedge(x, y, z, c)
%  
%  Input: 
%      x, y, z - the matrices for a surface
%      c - color of the curve. Default values is 'k'
%

if nargin==3
    c='k';
end
tf = ishold;
d=size(x); 
if length(d)~=2
    error('Input is not the matrices for a surface!');
end
[m,n]=size(x); 
r=x(:,1); s=y(:,1); t=z(:,1); plot3(r,s,t,c); hold on;
r=x(:,n); s=y(:,n); t=z(:,n); plot3(r,s,t,c);
r=x(1,:); s=y(1,:); t=z(1,:); plot3(r,s,t,c);
r=x(m,:); s=y(m,:); t=z(m,:); plot3(r,s,t,c);
if tf
    hold on;
else
    hold off;
end