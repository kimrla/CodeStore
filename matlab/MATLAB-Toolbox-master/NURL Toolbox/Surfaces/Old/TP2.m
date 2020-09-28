
% Get the parametric intersection curves
pts1=[x(1,:); x(2,:); zeros(1, length(dt))];
pts2=[x(3,:); x(4,:); zeros(1, length(dt))];
crv1=nulpts2crv(pts1, order);
crv2=nrlmake(pts2, crv1.knots, [0, 1], order);

% Insert more parametric points
n=length(dt);
dp=zeros(1, n-1); 
for i=1:n-1
    dp(i)=norm(pts1(:,i+1)-pts1(:,i));
end
ap=sum(dp)/(n-1); ip=find(dp>ap);
t=[]; 
for i=1:length(ip)
    ti=linspace(crv1.knots(ip(i)), crv1.knots(ip(i)+1), 3+round(dp(ip(i))/ap));
    ti([1,end])=[]; 
    t=[t, ti]; 
end
ps1=nrleval(crv1, t);
ps2=nrleval(crv2, t);
n=length(t);
xn=[]; pts1=[]; pts2=[]; dn=[];
% Get intersection points of two surfaces
for i=1:n
    % Extract a curve from a surface
    crv=nrlsrf2crv(srf2, [], ps2(2,i)); 

    % Get intersection points of the curve with the surface
    [xi, pnts1i, pnts2i, dti]=srfcrvintersects(srf1, crv, [ps1(1:2,i); ps2(1,i)]);
    xi(4,:)=ps2(2,i);
    
    % Assemble the results
    xn=[xn, xi]; dn=[dn; dti];
    pts1=[pts1, pnts1i];
    pts2=[pts2, pnts2i];
end
dc=dn<tol*1e-6;
xn=xn(:,dc); 
pts1=pts1(:,dc); 
pts2=pts2(:,dc); 
dn=dn(dc);

% Remove duplicated intersections
x=[x, xn]; dt=[dt; dn];
pnts1=[pnts1, pts1];
pnts2=[pnts2, pts2];
[x, pnts1, pnts2, dt]=optintersects(x, pnts1, pnts2, dt);

% Get the parametric intersection curves
pts1=[x(1,:); x(2,:); zeros(1, length(dt))];
pts2=[x(3,:); x(4,:); zeros(1, length(dt))];
crv1=nulpts2crv(pts1, order);
crv2=nrlmake(pts2, crv1.knots, [0, 1], order);


%% Least spuare approximation
% order=4; 
% [xx, i1]=min(x); [xx, i2]=min(xx); 
% ii=[1, 2, 3, 4]; jj=i1(i2); ii(jj)=[];
% [index, iu]=spanindex(x(jj,:), x(jj,:), order-1);
% [n, ~]=size(iu); y=x;
% for i=1:n
%     u=x(jj, iu(i,:));
%     A(1,1)=4; A(1,2)=sum(u);
%     A(2,1)=A(1,2);
%     A(2,2)=sum(u.^2);
%     for j=index(i)+1:index(i+1)
%         for k=1:3
%             v=x(ii(k), iu(i,:));
%             B(1,1) = sum(v);
%             B(2,1) = sum(u.*v);
%             a=A\B;
%             y(ii(k), j) = a(2)*x(jj, j)+a(1);
%         end
%     end
% end
% x=y;


