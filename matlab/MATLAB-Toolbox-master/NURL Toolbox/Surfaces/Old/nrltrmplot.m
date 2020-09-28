function nrltrmplot(crv1, srf1, crv, nuv, nt, pt)

if nargin==5
    pt=1;
end

ut1=linspace(0, 1, nt);
[pts1, jac1]=nrldeval(crv1, ut1);

tt1={linspace(0, 1, nuv(1)), linspace(0, 1, nuv(2))};
[v1, u1]=meshgrid(tt1{2}, tt1{1}); 
u1=u1(:); v1=v1(:); 
dm=DistanceMatrix([u1, v1], pts1(1:2,:)');
[~, I]=min(dm,[],2);
dr=[u1'; v1'; zeros(size(u1'))]-pts1(:,I);
dn=cross(dr, jac1(:,I));
p=dn(3,:)>0;
p=reshape(p, nuv(1), nuv(2));
p1=nrleval(srf1, tt1);
x=squeeze(p1(1,:,:)); 
y=squeeze(p1(2,:,:)); 
z=squeeze(p1(3,:,:)); 
if pt==1
    x(~p)=nan;
    y(~p)=nan;
    z(~p)=nan;
else
    x(p)=nan;
    y(p)=nan;
    z(p)=nan;
end

tf = ishold;
surf(x, y, z);
hold on;
colormap summer;
shading interp;
pnts=nrleval(crv, linspace(0,1,nt)); 
plot3(pnts(1,:), pnts(2,:), pnts(3,:), 'LineWidth', 2);
view(3);
if tf
    hold on;
else
    hold off;
end




