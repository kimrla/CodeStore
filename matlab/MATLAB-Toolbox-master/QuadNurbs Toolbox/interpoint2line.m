function [pts, x, d]=interpoint2line(L, P, x)

% interpoint2line: Get the nearest point on the line to the point
% 
% Calling Sequences:
% 
%     [pts,x, d]=interpoint2line(L, P, x)
% 
% INPUTS:
%
%      L=[P1x, P1y, P1z
%           P2x, P2y, P2z]
%           Two vertexes of a line:
% 
%      x - initial guess of parametric point on the two line. 
%           The default value is 0.5;
% 
% OUTPUT:
% 
%     pts - The nearest point on the line.
% 
%      x -  Parametric point on the line that is close to the point P. 
% 
%     d - Distance of the two points
% 

if nargin==2
    x=0.5;
end

pts=(1-x(1))*L(1,:)+x(1)*L(2,:);
jac=L(2,:)-L(1,:);
F(1,1)=dot(pts-P, jac); 
dF(1,1)=dot(jac, jac); 
x=x-dF\F; 
pts=(1-x(1))*L(1,:)+x(1)*L(2,:);
if nargout>2
    d=norm(P-pts);
end

%% demo
% % Create a line and a point
% L=[0,0,0; 1,0.5,0];
% P=[0.5,1,0];
% 
% % Get the nearest point on the line to the point
% [pts, x, d]=interpoint2line(L, P);
% 
% % Plot the results
% figure; hold on;
% plot3(L(:,1), L(:,2), L(:,3));
% plot3(P(1), P(2), P(3), 'ro');
% plot3(pts(1), pts(2), pts(3), 'r*');
% axis equal;




