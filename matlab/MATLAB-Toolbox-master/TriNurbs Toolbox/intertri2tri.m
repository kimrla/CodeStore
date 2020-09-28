function [pts1, pts2, x]=intertri2tri(T1, T2)

% interline2tri: Get the interections of two triangles
% 
% Calling Sequences:
% 
%     [pts1, pts2, x]=intertri2tri(T1, T2)
% 
% INPUTS:
%
%      T1, T2=[P1x, P1y, P1z
%                   P2x, P2y, P2z
%                   P3x, P3y, P3z]
%           Three vertexes of a triangle.
% 
% OUTPUT:
% 
%     pts1 - Points of intersections on the triangle.
% 
%     pts2 - Points of intersections on the line.
% 
%     x =[u1, v1, u2, v2], where [u, v] are parametric points of 
%          intersections for the triangles.
% 

[t1, tol1]=islinetri(T1);
[t2, tol2]=islinetri(T2);
if t1 || t2
    error('One of the two triangles is a straight line!');
end

t=0;
v1=cross(T1(1,:)-T1(2,:), T1(1,:)-T1(3,:));
v2=cross(T2(1,:)-T2(2,:), T2(1,:)-T2(3,:));
vt=norm(cross(v1, v2));
tol=max([tol1, tol2]); 
if vt<tol*1e-10    
    [pts1(1,:), x(1:2,1), d1]=interpoint2tri(T1, T2(1,:)); % Vertex 1 of triangle 2
    x(3,1)=1; x(4,1)=0; pts2(1,:)=T2(1,:); 
    [pts1(2,:), x(1:2,2), d2]=interpoint2tri(T1, T2(2,:)); % Vertex 2 of triangle 2
    x(3,2)=0; x(4,2)=1; pts2(2,:)=T2(2,:); 
    [pts1(3,:), x(1:2,3), d3]=interpoint2tri(T1, T2(3,:)); % Vertex 3 of triangle 2
    x(3,3)=0; x(4,3)=0; pts2(3,:)=T2(3,:); 
    if d1<tol*1e-10 && d2<tol*1e-10 && d3<tol*1e-10
        [pts2(4,:), x(3:4,4)]=interpoint2tri(T2, T1(1,:)); % Vertex 1 of triangle 1
        x(1,4)=1; x(2,4)=0; pts1(4,:)=T1(1,:); 
        [pts2(5,:), x(3:4,5)]=interpoint2tri(T2, T1(2,:)); % Vertex 2 of triangle 1
        x(1,5)=0; x(2,5)=1; pts1(5,:)=T1(2,:); 
        [pts2(6,:), x(3:4,6)]=interpoint2tri(T2, T1(3,:)); % Vertex 3 of triangle 1
        x(1,6)=0; x(2,6)=0; pts1(6,:)=T1(3,:); 
        t=6;
    end
end
[tpts1, tpts2, tx, dk]=interline2tri(T1, T2([1,2],:)); % Edge 1 of triangle 2
if dk<tol*1e-10
    tk=size(tx, 2); 
    x(1:2,t+1:t+tk)=tx(1:2,:); 
    x(3,t+1:t+tk)=1-tx(3,:); x(4,t+1:t+tk)=tx(3,:); 
    pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
    t=t+tk; 
end
[tpts1, tpts2, tx, dk]=interline2tri(T1, T2([3,1],:)); % Edge 2 of triangle 2
if dk<tol*1e-10
    tk=size(tx, 2); 
    x(1:2,t+1:t+tk)=tx(1:2,:); 
    x(3,t+1:t+tk)=tx(3,:); x(4,t+1:t+tk)=0; 
    pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
    t=t+tk; 
end
[tpts1, tpts2, tx, dk]=interline2tri(T1, T2([3,2],:)); % Edge 3 of triangle 2
if dk<tol*1e-10
    tk=size(tx, 2); 
    x(1:2,t+1:t+tk)=tx(1:2,:); 
    x(3,t+1:t+tk)=0; x(4,t+1:t+tk)=tx(3,:); 
    pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
    t=t+tk; 
end
[tpts1, tpts2, tx, dk]=interline2tri(T2, T1([1,2],:)); % Edge 1 of triangle 1
if dk<tol*1e-10
    tk=size(tx, 2); 
    x(3:4,t+1:t+tk)=tx(1:2,:); 
    x(1,t+1:t+tk)=1-tx(3,:); x(2,t+1:t+tk)=tx(3,:); 
    pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
    t=t+tk; 
end
[tpts1, tpts2, tx, dk]=interline2tri(T2, T1([3,1],:)); % Edge 2 of triangle 1
if dk<tol*1e-10
    tk=size(tx, 2); 
    x(3:4,t+1:t+tk)=tx(1:2,:); 
    x(1,t+1:t+tk)=tx(3,:); x(2,t+1:t+tk)=0; 
    pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
    t=t+tk; 
end
[tpts1, tpts2, tx, dk]=interline2tri(T2, T1([3,2],:)); % Edge 3 of triangle 1
if dk<tol*1e-10
    tk=size(tx, 2); 
    x(3:4,t+1:t+tk)=tx(1:2,:); 
    x(1,t+1:t+tk)=0; x(2,t+1:t+tk)=tx(3,:); 
    pts1(t+1:t+tk,:)=tpts1(:,:); pts2(t+1:t+tk,:)=tpts2(:,:); 
end
p=true(1,size(x,2)); 
p=p & x(1,:)>-eps & x(1,:)<1+eps; 
p=p & x(2,:)>-eps & x(2,:)<1+eps; 
p=p & x(1,:)+x(2,:)>-eps & x(1,:)+x(2,:)<1+eps; 
p=p & x(3,:)>-eps & x(3,:)<1+eps; 
p=p & x(4,:)>-eps & x(4,:)<1+eps; 
p=p & x(3,:)+x(4,:)>-eps & x(3,:)+x(4,:)<1+eps; 
x=x(:,p); 
pts1=pts1(p,:); 
pts2=pts2(p,:); 
[x, p]=RemDuplicate(x');
x=x';
pts1=pts1(p,:); 
pts2=pts2(p,:); 


%% demo 1 - two triangles are interected normally
% % Create a triangle
% T1=[0,0,0; 1,0,0; 1,1,1];
% T2=[0.5,1,0;0.5,0.2,0;0.6,0.1,0.8];
% 
% % Get the interections of two triangles
% [pts1, pts2, x]=intertri2tri(T1, T2);
% % pts1=tripoint(T1(1,:), T1(2,:), T1(3,:), x(1,:), x(2,:));
% % pts2=tripoint(T2(1,:), T2(2,:), T2(3,:), x(3,:), x(4,:));
% 
% % Plot the triangle
% n=15;
% tri=tridelaunay(n);
% [u, v]=TrigNodeVect(n);
% Pnts1=tripoint(T1(1,:), T1(2,:), T1(3,:), u, v);
% Pnts2=tripoint(T2(1,:), T2(2,:), T2(3,:), u, v);
% figure; hold on;
% trisurf(tri,Pnts1(:,1),Pnts1(:,2),Pnts1(:,3));
% trisurf(tri,Pnts2(:,1),Pnts2(:,2),Pnts2(:,3));
% plot3(pts1(:,1), pts1(:,2), pts1(:,3), 'ro');
% plot3(pts2(:,1), pts2(:,2), pts2(:,3), 'r*');
% view(3); axis equal;


%% demo 2 - One triangle is on the other one
% % Create a triangle
% T1=[0,0,0; 1,0,0; 1,1,1];
% P1=tripoint(T1(1,:), T1(2,:), T1(3,:), 0.1, 1.0);
% P2=tripoint(T1(1,:), T1(2,:), T1(3,:), 0.6, 0.3);
% P3=tripoint(T1(1,:), T1(2,:), T1(3,:), 0.1, 0.6);
% T2=[P1; P2; P3];
% 
% % Get the interections of two triangles
% [pts1, pts2, x]=intertri2tri(T1, T2);
% % pts1=tripoint(T1(1,:), T1(2,:), T1(3,:), x(1,:), x(2,:));
% % pts2=tripoint(T2(1,:), T2(2,:), T2(3,:), x(3,:), x(4,:));
% 
% % Plot the triangle
% n=15;
% tri=tridelaunay(n);
% [u, v]=TrigNodeVect(n);
% Pnts1=tripoint(T1(1,:), T1(2,:), T1(3,:), u, v);
% Pnts2=tripoint(T2(1,:), T2(2,:), T2(3,:), u, v);
% figure; hold on;
% trisurf(tri,Pnts1(:,1),Pnts1(:,2),Pnts1(:,3));
% trisurf(tri,Pnts2(:,1),Pnts2(:,2),Pnts2(:,3));
% plot3(pts1(:,1), pts1(:,2), pts1(:,3), 'ro');
% plot3(pts2(:,1), pts2(:,2), pts2(:,3), 'r*');
% view(3); axis equal;




