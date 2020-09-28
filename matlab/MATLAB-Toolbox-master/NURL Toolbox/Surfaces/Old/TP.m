clear; clc;

R=4;  
s1=0; s2=2*pi; t1=0; t2=pi; 
center=[5, 5, 4]; 
srf1 = nrlsphere(R, center, s1, s2, t1, t2); 
srf2 = nrltestsrf; 
% srf2=nrl4surf([10.0 0.0 5.3], [0.0 0.0 3.3], [0.0 10.0 5.3], [10.0 10.0 3.3]);
figure; hold on; 
nrlplot(srf1, [100, 100]); 
nrlplot(srf2, [100, 100]); 
axis equal; view(3); 
shading interp; 

% First and second derivatives
m1=16; n1=16; m2=16; n2=17; 
s1=linspace(0, 1, m1); t1=linspace(0, 1, n1); 
s2=linspace(0, 1, m2); t2=linspace(0, 1, n2); 
tt1={s1, t1}; tt2={s2, t2}; 
[pnts1, jac1, hess1]=nrldeval(srf1, tt1); 
[pnts2, jac2, hess2]=nrldeval(srf2, tt2); 
[H1, K1]=srfcurvatures(jac1, hess1);
dt1=real(sqrt(H1.^2-K1));
k1=H1+dt1;
[H2, K2]=srfcurvatures(jac2, hess2);
dt2=real(sqrt(H2.^2-K2));
k2=H2+dt2;
[~, I]=min([-min(k1), -min(k2)]);

% Get tolerance of the surfaces
crvs1 = nrlextract(srf1); 
crvs2 = nrlextract(srf2); 
ln1=max([nrlmeasure(crvs1(1)), nrlmeasure(crvs1(2)), nrlmeasure(crvs1(3)), nrlmeasure(crvs1(4))]);
ln2=max([nrlmeasure(crvs2(1)), nrlmeasure(crvs2(2)), nrlmeasure(crvs2(3)), nrlmeasure(crvs2(4))]);
tol=max([ln1, ln2])/min([m1, n1, m2, n2]); 

% Get approximate intersection points
dm=DistanceMatrix(pnts1(:,:)', pnts2(:,:)'); 
pp=dm<tol;
qq=zeros(m1*n1, 1);
for i=1:m1*n1
    qi=find(pp(i,:));
    if ~isempty(qi)
        di=DistanceMatrix(pnts1(:,i)', pnts2(:,qi)'); 
        [~, id]=min(di);
        qq(i)=qi(id);
    end
end
id=qq~=0;
p=find(id); q=qq(id);

% Sort the points
ps=zeros(size(p)); qs=ps; m=length(ps);
if I==1
    dm=DistanceMatrix(pnts1(:,p)', pnts1(:,p)'); 
else
    dm=DistanceMatrix(pnts2(:,q)', pnts2(:,q)'); 
end
ps(1)=p(1); qs(1)=q(1);
id=true(m,1);
ip0=1;
for i=2:m
    id(ip0)=false;
    dd=dm(ip0, id);
    [~, ip] = min(dd);
    pp=find(id); 
    ip=pp(ip);
    ps(i)=p(ip); qs(i)=q(ip);
    ip0=ip;
end

% Remove duplicated points
figure; hold on; 
nrlplot(srf1, [100, 100]); 
shading interp; 
plot3(pnts1(1,ps), pnts1(2,ps), pnts1(3,ps));
axis equal; view(3); 

figure; hold on; 
nrlplot(srf2, [100, 100]); 
shading interp; 
axis equal; view(3); 
plot3(pnts2(1,qs), pnts2(2,qs), pnts2(3,qs));

% Get parametric points
[v1, u1]=meshgrid(t1, s1); u1=u1(:); v1=v1(:); 
[v2, u2]=meshgrid(t2, s2); u2=u2(:); v2=v2(:); 
u1=u1(ps); v1=v1(ps); 
u2=u2(qs); v2=v2(qs); 

% % Get intersection of two surfaces
% X=zeros(m, 4);
% pnts1=zeros(3,m); pnts2=pnts1;
% for k=1:m
%     x=[u1(k), v1(k), u2(k), v2(k)]';
%     [X(k,:), pnts1(:,k), pnts2(:,k)]=srfintersect(srf1, srf2, x);
% end

% figure; hold on; 
% nrlplot(srf1, [100, 100]); 
% shading interp; 
% plot3(pnts1(1,:), pnts1(2,:), pnts1(3,:), 'LineWidth',2);
% axis equal; view(3); 
% 
% figure; hold on; 
% nrlplot(srf2, [100, 100]); 
% shading interp; 
% plot3(pnts2(1,:), pnts2(2,:), pnts2(3,:), 'LineWidth',2);
% axis equal; view(3); 


%% Mesh generation
% srf = nrltestsrf; 
% figure; hold on; 
% nrlplot(srf, [100, 100]); 
% axis equal; view(3); 
% shading interp; 
% 
% % First and second derivatives
% m=16; n=17; 
% s=linspace(0, 1, m); t=linspace(0, 1, n); 
% tt={s, t}; 
% [pnts, jac, hess]=nrldeval(srf, tt); 
% 
% % First and second quadratic form
% [E, F, G, g, nor]=quadratic1(jac);
% [L, M, N]=quadratic2(hess, nor);
% 
% % Surface standard expansion at a point
% k=1; len=3; 
% r1=vecnorm(jac{1}(:,k)); nv=nor(:,k); r2=cross(nv, r1);
% quiver3(pnts(1,k), pnts(2,k), pnts(3,k), len*r1(1), len*r1(2), len*r1(3));
% quiver3(pnts(1,k), pnts(2,k), pnts(3,k), len*r2(1), len*r2(2), len*r2(3));
% quiver3(pnts(1,k), pnts(2,k), pnts(3,k), len*nv(1), len*nv(2), len*nv(3));


