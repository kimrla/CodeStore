function [pts1, pts2, x, d, inter]=interline2quad(Q, L, x)

% interline2quad: Get intersections of a line with a quadrangle
% 
% Calling Sequences:
% 
%     [pts1, pts2, x, d, inter]=interline2quad(Q, L, x)
% 
% INPUTS:
%
%      Q=[P1x, P1y, P1z
%           P2x, P2y, P2z
%           P3x, P3y, P3z
%           P4x, P4y, P4z]
%           Three vertexes of a quadrangle.
%
%      L=[P1x, P1y, P1z
%           P2x, P2y, P2z]
%           Two vertexes of a line:
% 
%      x =[u; v; s], where [u; v] are initial guess parametric points of 
%           intersection for the quadrangle, and s is that of the line.
%           The default values are [0.5; 0.5; 0.5];
% 
% OUTPUT:
% 
%     pts1 - Points of intersections on the quadrangle.
% 
%     pts2 - Points of intersections on the line.
% 
%     x =[u, v, s], where [u, v] are parametric points of 
%          intersection for the quadrangle, and s is that of the line.
% 
%     d - Distance of the two points
%
%     inter - A value that indicates whether the line is parallel
%              (inter==0) or intersected with the triangle (inter>0). 
%
% Discription:
%
%     The distance of the middle point of the line to the quadrangle
%     is returned if the line is parallel to the quadrangle.
% 

if nargin==2
    x=[0.5; 0.5; 0.5];
end
[t, tol]=islinetri(Q);
if t
    error('The triangle is a straight line!');
end
x=x(:);
for j=1:3
    [pts1, jac1]=quadpoint(Q, x(1), x(2));
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
        pts1=quadpoint(Q, x(1), x(2)); 
        pts2=(1-x(3))*L(1,:)+x(3)*L(2,:); 
        d=norm(pts2-pts1);

        % The case that the line is one the triangle
        if d<tol*1e-10
            % Two edges of the line
            x(3,1)=0; pts2(1,:)=L(1,:); d=zeros(1,2);
            [pts1(1,:), x(1:2,1), d(1)]=interpoint2quad(Q, pts2(1,:)); 
            x(3,2)=1; pts2(2,:)=L(2,:); 
            [pts1(2,:), x(1:2,2), d(2)]=interpoint2quad(Q, pts2(2,:)); 
            
            % Edge 1 of the quadrangle
            [tpts1, tpts2, tx, tdk]=interline2line(Q([1,2],:), L); 
            t=2; 
            tk=size(tx, 2); 
            pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
            x(1,t+1:t+tk)=1-tx(1); x(2,t+1:t+tk)=tx(1); x(3,t+1:t+tk)=tx(2); 
            d(t+1:t+tk)=tdk;
            t=t+tk; 
            
            % Edge 2 of the quadrangle
            [tpts1, tpts2, tx, tdk]=interline2line(Q([2,3],:), L); 
            tk=size(tx, 2); 
            pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
            x(1,t+1:t+tk)=1-tx(1); x(2,t+1:t+tk)=tx(1); x(3,t+1:t+tk)=tx(2); 
            d(t+1:t+tk)=tdk;
            t=t+tk; 
            
            % Edge 3 of the quadrangle
            [tpts1, tpts2, tx, tdk]=interline2line(Q([3,4],:), L); 
            tk=size(tx, 2); 
            pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
            x(1,t+1:t+tk)=1-tx(1); x(2,t+1:t+tk)=tx(1); x(3,t+1:t+tk)=tx(2); 
            d(t+1:t+tk)=tdk;
            t=t+tk; 
            
            % Edge 4 of the quadrangle
            [tpts1, tpts2, tx, tdk]=interline2line(Q([4,1],:), L); 
            tk=size(tx, 2); 
            pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
            x(1,t+1:t+tk)=1-tx(1); x(2,t+1:t+tk)=tx(1); x(3,t+1:t+tk)=tx(2); 
            d(t+1:t+tk)=tdk;
            break;
        end
    else
        F(3,1)=-dot(dr, jac2);     
        dF(1,3)=-dot(jac1{1}, jac2); 
        dF(3,1)=dF(1,3);     
        dF(2,3)=-dot(jac1{2}, jac2); 
        dF(3,2)=dF(2,3); 
        dF(3,3)=dot(jac2, jac2); 
        x=x-dF\F; 
        pts1=quadpoint(Q, x(1), x(2));
        pts2=(1-x(3))*L(1,:)+x(3)*L(2,:);
        d=norm(pts2-pts1);
    end
end
p=x(3,:)>-eps & x(3,:)<1+eps;
p=p & d<tol*1e-3;
p=p & x(1,:)>-eps & x(1,:)<1+eps; 
p=p & x(2,:)>-eps & x(2,:)<1+eps; 
x=x(:,p); 
pts1=pts1(p,:); 
pts2=pts2(p,:); 
[x, p]=RemDuplicate(x');
x=x';
pts1=pts1(p,:); 
pts2=pts2(p,:); 


%% demo
% % The points
% L=[0.0,0.1,0.8; 1,1,0.8];
% Q=[0,0,1; 1,0,0; 1,1.0,1; 0,1,0];
% 
% % Get the inter section of the line to the quandrangle
% [pe, pq, x, d]=interline2quad(Q, L, [0,0,0]);
% 
% % Test of quadpoint
% m=5; n=6;
% s=linspace(0,1,m);
% t=linspace(0,1,n);
% [t, s]=meshgrid(t, s);
% [Pt, Jac, Hess]=quadpoint(Q, s, t);
% Px=reshape(Pt(:,1), m, n);
% Py=reshape(Pt(:,2), m, n);
% Pz=reshape(Pt(:,3), m, n);
% 
% % Plot results
% figure; surf(Px, Py, Pz);
% hold on;
% plot3(L(:,1), L(:,2), L(:,3), 'b', 'LineWidth', 2);
% plot3(pe(:,1), pe(:,2), pe(:,3), 'm*');
% plot3(pq(:,1), pq(:,2), pq(:,3), 'ro');
% axis equal;




