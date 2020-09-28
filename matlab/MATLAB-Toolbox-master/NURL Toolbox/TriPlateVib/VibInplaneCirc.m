% Test triangular element
clear; clc;

% Young's modulus (E), Poisson's ratio (nu), density (rho) and radious (R)
% Order of the NURLS basis (order), the number of nodes (nodes)
E=71e9; nu=0.33; rho=2700; R=0.5;
order=3*[1, 1]; nodes=5*[1, 1];

% Node numbers of elements
num1=nodes(2); num2=nodes(1); 
[Num, BN, TN]=NodeNumTri(num1, num2);
TNe=num1*(num2-1)+1;
ss=1:TNe; tt=zeros(1, 2*TNe);
s1=2*ss-1; s2=2*ss;
tt(s1)=ss(:); tt(s2)=TNe+ss(:);

% Get element matrices
elemat=cell(4,1);
figure; hold on;
for j=1:4
    % Creat a triangular nurl patch
    tsrf=SectorGeom(R, j, order, nodes);
    elemat{j}.tsrf=tsrf;

    % Stiffness and mass matrices
    [Kj, Mj, pnts, elemat{j}.G]=StiffMassInpTri(tsrf, E, nu, rho);
    elemat{j}.pts=tsrf.faces.coefs(1:3,1:TNe);
    elemat{j}.Ke=Kj(tt, tt); 
    elemat{j}.Me=Mj(tt, tt); 
    
    mn=size(pnts); m=mn(2); n=mn(3);
    elemat{j}.X=squeeze(pnts(1,:,:));
    elemat{j}.Y=squeeze(pnts(2,:,:));
    nrlctrlplot(tsrf.faces);
end

% Plot node numbers
figure; hold on;
for j=1:4
    surf(elemat{j}.X, elemat{j}.Y, 0*elemat{j}.X);
    plot(elemat{j}.X(:,[1,end]), elemat{j}.Y(:,[1,end]), 'r');
    plot(elemat{j}.X([1,end],:)', elemat{j}.Y([1,end],:)', 'r');
    for i=1:TNe
        text(elemat{j}.pts(1, i), elemat{j}.pts(2, i), num2str(Num(j, i)));
    end
end
axis equal;
shading interp;

% Assemble elements
Ke=zeros(2*TN); Me=Ke;
for j=1:4
    Ke=AssembleInP(Ke, elemat{j}.Ke, Num(j,:)); 
    Me=AssembleInP(Me, elemat{j}.Me, Num(j,:)); 
end

% Apply boundary conditions
bc=true(1, 2*TN); 
bc(2*BN-1)=false;
bc(2*BN)=false;
Ke=Ke(bc,bc);
Me=Me(bc,bc);

% Solve Eigenvalues (SE): 1 - all, 2,3- first a few
Ns=20; SE=3;
[d, V, NT]=SolveEig(Ke, Me, Ns, SE);
[sd, I]=sort(d); 
Gs=E/(2*(1+nu));
Frq=real(sqrt(sd'))/(2*pi);
sd=real(sqrt(sd*rho/Gs));

% Draw mode shapes
k=4; dd=0.05;
Vt=zeros(2*TN, 1);
Vt(bc)=real(V(:, I(k)));
mv=0;
for j=1:4
    s1=2*Num(j,:)-1;
    s2=2*Num(j,:);
    Uv=Vt(s1); Vv=Vt(s2);
    Uv=elemat{j}.G*Uv; 
    Vv=elemat{j}.G*Vv; 
    Uv=reshape(Uv, m, n); 
    Vv=reshape(Vv, m, n); 
    elemat{j}.Uv=Uv;
    elemat{j}.Vv=Vv;
    mj=max(max(abs([Uv,Vv])));
    if mj>mv
        mv=mj;
    end
end
figure; hold on;
for j=1:4
    surf(elemat{j}.X,elemat{j}.Y,elemat{j}.Uv);
end
title('U'); view(3);
figure; hold on;
for j=1:4
    surf(elemat{j}.X,elemat{j}.Y,elemat{j}.Vv);
end
title('V'); view(3);

% Draw in-plane deformation
figure; hold on;
for j=1:4
    Xd=elemat{j}.X+dd*elemat{j}.Uv/mv;
    Yd=elemat{j}.Y+dd*elemat{j}.Vv/mv;
    plot(Xd, Yd);
    plot(Xd', Yd');
%     plot(elemat{j}.X(:,1), elemat{j}.Y(:,1), 'r-.');
end
axis equal;
axis off;



