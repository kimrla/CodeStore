function [pts1, pts2, x, d, inter]=interline2tri(T, L, x)

% interline2tri: Get intersections of a line with a triangle
% 
% Calling Sequences:
% 
%     [pts1, pts2, x, d, inter]=interline2tri(T, L)
% 
% INPUTS:
%
%      T=[P1x, P1y, P1z
%           P2x, P2y, P2z
%           P3x, P3y, P3z]
%           Three vertexes of a triangle.
%
%      L=[P1x, P1y, P1z
%           P2x, P2y, P2z]
%           Two vertexes of a line:
% 
%      x =[u; v; s], where [u; v] are initial guess parametric points of 
%           intersection for the triangle, and s is that of the line.
%           The default values are [0.5; 0.5; 0.5];
% 
% OUTPUT:
% 
%     pts1 - Points of intersections on the triangle.
% 
%     pts2 - Points of intersections on the line.
% 
%     x =[u, v, s], where [u, v] are parametric points of 
%          intersection for the triangle, and s is that of the line.
% 
%     d - Distance of the line to the triangle.
%
%     inter - A value that indicates whether the line is parallel
%              (inter==0) or intersected with the triangle (inter>0). 
%
% Discription:
%
%     The distance of the middle point of the line to the triangle
%     is returned if the line is parallel to the triangle.
% 

if nargin==2
    x=[0.5; 0.5; 0.5];
end
[t, tol]=islinetri(T);
if t
    error('The triangle is a straight line!');
end
x=x(:);
[pts1, jac1]=tripoint(T(1,:), T(2,:), T(3,:), x(1), x(2));
pts2=(1-x(3))*L(1,:)+x(3)*L(2,:);
jac2=L(2,:)-L(1,:);
inter=abs(dot(cross(jac1{1}, jac1{2}), jac2));
dr=pts1-pts2; 
F(1,1)=dot(dr, jac1{1}); 
F(2,1)=dot(dr, jac1{2}); 
dF(1,1)=dot(jac1{1}, jac1{1}); 
dF(1,2)=dot(jac1{1}, jac1{2}); 
dF(2,1)=dF(1,2);
dF(2,2)=dot(jac1{2}, jac1{2}); 
if inter<tol*1e-10
    x(1:2)=x(1:2)-dF\F; 
    pts1=tripoint(T(1,:), T(2,:), T(3,:), x(1), x(2)); 
    pts2=(1-x(3))*L(1,:)+x(3)*L(2,:); 
    d=norm(pts2-pts1); 
    
    % The case that the line is one the triangle
    if d<tol*1e-10
        x(3,1)=0; pts2(1,:)=L(1,:); d=zeros(1,2);
        [pts1(1,:), x(1:2,1), d(1)]=interpoint2tri(T, pts2(1,:)); 
        x(3,2)=1; pts2(2,:)=L(2,:); 
        [pts1(2,:), x(1:2,2), d(2)]=interpoint2tri(T, pts2(2,:)); 
        
        % Edge 1 of the triangle
        [tpts1, tpts2, tx, tdk]=interline2line(T([1,2],:), L); 
        t=2; 
        tk=size(tx, 2); 
        pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
        x(1,t+1:t+tk)=1-tx(1); x(2,t+1:t+tk)=tx(1); x(3,t+1:t+tk)=tx(2); 
        d(t+1:t+tk)=tdk;
        t=t+tk; 
            
        % Edge 2 of the triangle
        [tpts1, tpts2, tx, tdk]=interline2line(T([3,1],:), L); 
        tk=size(tx, 2); 
        pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
        x(1,t+1:t+tk)=tx(1); x(2,t+1:t+tk)=0; x(3,t+1:t+tk)=tx(2); 
        d(t+1:t+tk)=tdk;
        t=t+tk; 

        % Edge 3 of the triangle
        [tpts1, tpts2, tx, tdk]=interline2line(T([3,2],:), L); 
        tk=size(tx, 2); 
        pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
        x(1,t+1:t+tk)=0; x(2,t+1:t+tk)=tx(1); x(3,t+1:t+tk)=tx(2); 
        d(t+1:t+tk)=tdk;
        
        p=true(1,size(x,2)); 
        p=p & d<tol*1e-10;
        p=p & x(3,:)>-eps & x(3,:)<1+eps; 
        p=p & x(1,:)>-eps & x(1,:)<1+eps; 
        p=p & x(2,:)>-eps & x(2,:)<1+eps; 
        p=p & x(1,:)+x(2,:)>-eps & x(1,:)+x(2,:)<1+eps; 
        x=x(:,p); 
        pts1=pts1(p,:); 
        pts2=pts2(p,:); 
        d=d(:,p); 
        [x, p]=RemDuplicate(x');
        x=x';
        pts1=pts1(p,:); 
        pts2=pts2(p,:); 
        d=d(:,p); 
    end
else
    F(3,1)=-dot(dr, jac2);     
    dF(1,3)=-dot(jac1{1}, jac2); 
    dF(3,1)=dF(1,3);     
    dF(2,3)=-dot(jac1{2}, jac2); 
    dF(3,2)=dF(2,3); 
    dF(3,3)=dot(jac2, jac2); 
    x=x-dF\F; 
    pts1=tripoint(T(1,:), T(2,:), T(3,:), x(1), x(2));
    pts2=(1-x(3))*L(1,:)+x(3)*L(2,:);
    d=norm(pts2-pts1);
end


%% demo 1 - the line is not on the triangle
% % Create a triangle
% T=[0,0,0; 1,0,0; 1,1,1];
% L=[0.5,0.5,0; 0.9,0.3,1];
% % L=[0,0.5,0; 1,0.5,0];
% 
% % Get the interection
% [pts1, pts2, x, d]=interline2tri(T, L);
% 
% % Plot the results
% n=10;
% tri=tridelaunay(n);
% [u, v]=TrigNodeVect(n);
% Ps=tripoint(T(1,:), T(2,:), T(3,:), u, v);
% figure; hold on;
%  trisurf(tri, Ps(:,1),Ps(:,2),Ps(:,3));
% plot3(L(:,1), L(:,2), L(:,3));
% plot3(pts1(1), pts1(2), pts1(3), 'ro');
% plot3(pts2(1), pts2(2), pts2(3), 'r*');
% axis equal; view(3); 


%% demo 2 - the line is on the triangle
% % Create a triangle
% T=[0,0,0; 1,0,0; 1,1,1];
% P1=tripoint(T(1,:), T(2,:), T(3,:), 0.6, 0.4);
% P2=tripoint(T(1,:), T(2,:), T(3,:), 0.4, 0.6);
% L=[P1; P2];
% 
% % Get the interection
% [pts1, pts2, x, d, inter]=interline2tri(T, L);
% 
% % Plot the triangle
% n=10;
% tri=tridelaunay(n);
% [u, v]=TrigNodeVect(n);
% Ps=tripoint(T(1,:), T(2,:), T(3,:), u, v);
% figure; hold on;
%  trisurf(tri, Ps(:,1),Ps(:,2),Ps(:,3));
% plot3(L(:,1), L(:,2), L(:,3), 'r', 'LineWidth', 2);
% plot3(pts1(:,1), pts1(:,2), pts1(:,3), 'ro');
% plot3(pts2(:,1), pts2(:,2), pts2(:,3), 'r*');
% axis equal; view(3); 



