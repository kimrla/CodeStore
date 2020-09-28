function  [x, pnts1, pnts2, dt, tol]=srfsintersct(srf1, srf2, it, pts1, pnts1, jac1, pts2, pnts2, jac2)

% Get intersection points of two surfaces
% 
% Calling Sequences:
% 
%     [x, pnts1, pnts2, dt]=srfsintersct(srf1, srf2) 
%     [x, pnts1, pnts2, dt]=srfsintersct(srf1, srf2, it) 
%     [x, pnts1, pnts2, dt]=srfsintersct(srf1, srf2, pts1, pnts1, jac1, pts2, pnts2, jac2) 
%     [x, pnts1, pnts2, dt]=srfsintersct(srf1, srf2, it, pts1, pnts1, jac1, pts2, pnts2, jac2) 
% 
% INPUTS:
%
%     srf1, srf2 - nurls surfaces 
%     it - the number of iterations 
%     pts1 - approximated parametric points of intersection for srf1 
%     pnts1 - approximated points of intersection for srf1 
%     jac1 - first order derivatives on pnts1  
%     pts2 - approximated parametric points of intersection for srf2 
%     pnts2 - approximated points of intersection for srf2 
%     jac2 - first order derivatives on pnts2 
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     pnts1 - points of intersections on surface 1
%     pnts2 - points of intersections on surface 2
%     dt - distance of the two nearest points
%     tol - tolerance of length
% 
%  See also:
%       srfsinterscts
% 

% Get nearest distances of two surfaces and the corresponding derivatives
if nargin<=3
    % Tolerance of length and sampling points on the surfaces
    [pnts1, tt1, tol1, jac1]=srftolerance(srf1);
    [pnts2, tt2, tol2, jac2]=srftolerance(srf2);
    tol=max([tol1, tol2]);
    [p, q]=apntsints(pnts1, pnts2, tol); 

    % Get intersection points and corresponding parametric points
    [v1, u1]=meshgrid(tt1{2}, tt1{1}); u1=u1(:); v1=v1(:); 
    [v2, u2]=meshgrid(tt2{2}, tt2{1}); u2=u2(:); v2=v2(:); 
    pts1=[u1(p)'; v1(p)']; 
    pts2=[u2(q)'; v2(q)']; 
    pnts1=pnts1(:,:); pnts1=pnts1(:,p); 
    jac1{1}=jac1{1}(:,:); jac1{1}=jac1{1}(:,p); 
    jac1{2}=jac1{2}(:,:); jac1{2}=jac1{2}(:,p); 
    pnts2=pnts2(:,:); pnts2=pnts2(:,q); 
    jac2{1}=jac2{1}(:,:); jac2{1}=jac2{1}(:,q); 
    jac2{2}=jac2{2}(:,:); jac2{2}=jac2{2}(:,q); 
    if nargin==2
        it=4;
    end
end

if nargin==8
    jac2=pnts2;
    pnts2=pts2;
    pts2=jac1;
    jac1=pnts1;
    pnts1=pts1;
    pts1=it;
end

% Get maximum angles of the tangent vectors of approximated intersection points
n=length(pts1(1,:)); ag=zeros(n,4);
for i=1:n
    for j=1:2
        v1=jac1{j}(:, i); v2=jac2{1}(:, i); 
        ag(i, 2*(j-1)+1)=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
        v1=jac1{j}(:, i); v2=jac2{2}(:, i); 
        ag(i, 2*(j-1)+2)=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
    end
end
q=ag>pi/2; 
ag(q)=pi-ag(q);
[~, I] = max(ag, [], 2);

% Get improved approximated intersection points
ii=[1, 1, 2, 2; 2, 2, 1, 1]; jj=[1, 2, 1, 2; 2, 1, 2, 1];
dt=zeros(n,1);
for j=1:it
    % Use nearest points of two straight lines
    d=mod(j, 2);
    if d==0
        d=2;
    end
    for i=1:n
        i1=ii(d, I(i)); i2=jj(d, I(i));
        p1=pnts1(:,i); dp1=jac1{i1}(:,i);
        p2=pnts2(:,i); dp2=jac2{i2}(:,i);
        x=[pts1(i1, i); pts2(i2, i)]; 
        F=[dot(p1-p2, dp1); -dot(p1-p2, dp2)];
        dF=[dot(dp1, dp1), -dot(dp1, dp2)
              -dot(dp1, dp2),  dot(dp2, dp2)];
        x=x-dF\F;
        p=x>1; x(p)=1; p=x<0; x(p)=0; 
        pts1(i1, i)=x(1); 
        pts2(i2, i)=x(2);
    end
    
    % Use nearest points of points to a plane
    [pnts1, jac1]=nrldeval(srf1, pts1); 
    [pnts2, jac2]=nrldeval(srf2, pts2); 
    for i=1:n
        p1=pnts1(:,i); p2=pnts2(:,i); 
        dp21=jac2{1}(:,i); 
        dp22=jac2{2}(:,i); 
        x=pts2(:,i); 
        F=[dot(p2-p1, dp21); dot(p2-p1, dp22)]; 
        dF=[dot(dp21, dp21), dot(dp21, dp22) 
               dot(dp21, dp22), dot(dp22, dp22)]; 
        x=x-dF\F; 
        pts2(:, i)=x(:); 
    end
    p=pts2>1; pts2(p)=1; p=pts2<0; pts2(p)=0; 
    
    % Use nearest points of points to a plane
    [pnts2, jac2]=nrldeval(srf2, pts2);
    for i=1:n
        p1=pnts1(:,i); p2=pnts2(:,i);  
        dp11=jac1{1}(:,i); 
        dp12=jac1{2}(:,i); 
        x=pts1(:,i); 
        F=[dot(p1-p2, dp11); dot(p1-p2, dp12)]; 
        dF=[dot(dp11, dp11), dot(dp11, dp12) 
               dot(dp11, dp12), dot(dp12, dp12)]; 
        x=x-dF\F;
        pts1(:, i)=x(:);
    end
    p=pts1>1; pts1(p)=1; p=pts1<0; pts1(p)=0; 
    [pnts1, jac1]=nrldeval(srf1, pts1);     
end
for i=1:n
    dt(i)=norm(pnts1(:,i)-pnts2(:,i));
end
x=[pts1; pts2];




