function [pts1, pts2, x, d, inter]=interline2line(L1, L2, x)

% interline2line: Get intersections of two lines
% 
% Calling Sequences:
% 
%     [pts1, pts2, x, d, inter]=interline2line(L1, L2, x)
% 
% INPUTS:
%
%      L1, L2=[P1x, P1y, P1z
%           P2x, P2y, P2z]
%           Two vertexes of a line:
% 
%      x =[u1, u2], where u1 and u2 are respectively 
%           initial guess parametric points of intersection 
%           the two lines. The default values are [0.5; 0.5];
% 
% OUTPUT:
% 
%     pts1 - Points of intersections on the triangle.
% 
%     pts2 - Points of intersections on the line.
% 
%      x =[u1, u2], where u1 and u2 are respectively 
%           parametric points of intersection the two lines. 
% 
%     d - Distance of the two points
%
%     inter - A value that indicates whether the line is parallel
%              (inter==0) or intersected with the triangle (inter>0). 
%
% Discription:
%
%     The distance and the interections of end points
%     are returned if the two lines are parallel.
% 

if nargin==2
    x=[0.5; 0.5];
end

x=x(:);
jac1=L1(2,:)-L1(1,:);
jac2=L2(2,:)-L2(1,:);
inter=norm(cross(jac1, jac2));
tol=max([norm(jac1), norm(jac2)]);
if inter<tol*1e-10
    pts1=zeros(4,3); pts2=pts1; x=zeros(2,4); 
    d=zeros(1,4);
    pts1(1,:)=L1(1,:); x(1,1)=0;
    [pts2(1,:), x(2,1), d(1)]=interpoint2line(L2, pts1(1,:));
    pts1(2,:)=L1(2,:); x(1,2)=1;
    [pts2(2,:), x(2,2), d(2)]=interpoint2line(L2, pts1(2,:));
    pts2(3,:)=L2(1,:); x(2,3)=0;
    [pts1(3,:), x(1,3), d(3)]=interpoint2line(L1, pts2(3,:));
    pts2(4,:)=L2(2,:); x(2,4)=1;
    [pts1(4,:), x(1,4), d(4)]=interpoint2line(L1, pts2(4,:));
    return;    
else
    pts1=(1-x(1))*L1(1,:)+x(1)*L1(2,:);
    pts2=(1-x(2))*L2(1,:)+x(2)*L2(2,:);
    dr=pts1-pts2; 
    F(1,1)=dot(dr, jac1);     
    F(2,1)=-dot(dr, jac2);     
    dF(1,1)=dot(jac1, jac1); 
    dF(1,2)=-dot(jac1, jac2); 
    dF(2,1)=dF(1,2);  
    dF(2,2)=dot(jac2, jac2); 
    x=x-dF\F; 
end
pts1=(1-x(1))*L1(1,:)+x(1)*L1(2,:);
pts2=(1-x(2))*L2(1,:)+x(2)*L2(2,:);
d=norm(pts2-pts1);


%% demo 1 - interected
% % Create two lines
% L1=[0,0.5,0; 1,0.5,0];
% L2=[0,0,0; 1,1,0];
% 
% % Get interections of the two lines
% [pts1, pts2, x, d, inter]=interline2line(L1, L2);
% 
% % Plot the results
% figure; hold on;
% plot3(L1(:,1), L1(:,2), L1(:,3));
% plot3(L2(:,1), L2(:,2), L2(:,3));
% plot3(pts1(1), pts1(2), pts1(3), 'ro');
% plot3(pts2(1), pts2(2), pts2(3), 'r*');
% view(2);

%% demo 2 - parallel
% % Create two lines
% L1=[0,0.5,0; 1,0.5,0];
% L2=[0.2,0,0; 0.8,0,0];
% 
% % Get interections of the two lines
% [pts1, pts2, x, d, inter]=interline2line(L1, L2);
% 
% % Plot the results
% Pts1=zeros(length(d), 3);
% Pts2=Pts1;
% for j=1:length(d)
%     Pts1(j,:)=(1-x(1,j))*L1(1,:)+x(1,j)*L1(2,:);
%     Pts2(j,:)=(1-x(2,j))*L2(1,:)+x(2,j)*L2(2,:);
% end
% figure; hold on;
% plot3(L1(:,1), L1(:,2), L1(:,3));
% plot3(L2(:,1), L2(:,2), L2(:,3));
% plot3(pts1(:,1), pts1(:,2), pts1(:,3), 'ro');
% plot3(Pts1(:,1), Pts1(:,2), Pts1(:,3), 'r*');
% plot3(pts2(:,1), pts2(:,2), pts2(:,3), 'ro');
% plot3(Pts2(:,1), Pts2(:,2), Pts2(:,3), 'r*');
% view(2);






