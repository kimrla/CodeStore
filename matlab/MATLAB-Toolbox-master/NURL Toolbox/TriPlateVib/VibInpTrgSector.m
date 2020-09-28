% Test triangular element
clear; clc;

% Young's modulus (E), Poisson's ratio (nu) and density (rho)
% Order of the NURLS basis (order), the number of nodes (nodes)
E=195e9; nu=0.3; rho=7722.7;
order=5*[1, 1]; nodes=10*[1, 1];

% Boundary conditions (bu, bv): 1-clamped, 0-free
bu=[0, 1, 0]; bv=[0, 0, 1];

% Creat a triangular nurls patch
crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
af=pi/2; a=cos(af); b=sin(af);
crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
crv3 = nrlcirc(1, [0, 0], 0, af);
srf = nrltrgcoons(crv1, crv2, crv3);
srf = nrldegelev(srf, order-srf.order);
iknts=nodes-srf.number;
srf=nrlkntins(srf, {iknts(1), iknts(2)});

% Transform the nurls triangle patch for subsequent manipulations
tsrf=nrl2trg(srf);
figure; hold on;
nrlctrlplot(tsrf.faces); view(2);
nrlaxisplot(tsrf.faces);

% The total number of nodes
srf=tsrf.faces;
num1=srf.number(1); 
num2=srf.number(2); 
TN=num1*(num2-1)+1;

% Get first derivatives and weights of derivatives in Cartesian coordinates
[G, Gx, Gy, S, T, C, pnts, Jacb]=nrlplanemattrg(tsrf);
mn=size(pnts); m=mn(2); n=mn(3);

% Stiffness and mass matrices
CJ=diag(C(:).*Jacb(:));
Cr=E/(1-nu^2); nu1=(1-nu)/2;
G11=Gx'*CJ*Gx;
G22=Gy'*CJ*Gy;
G12=Gx'*CJ*Gy;
Ke=Cr*[G11+nu1*G22, nu*G12+nu1*G12';
           nu*G12'+nu1*G12, G22+nu1*G11];
Gm=G'*CJ*G; 
Zm=zeros(size(Gm)); 
Me=rho*[Gm, Zm; Zm, Gm]; 

% Apply boundary conditions
bcu=true(1, TN); bcv=bcu;
bn={1:num1, [num1:num1:TN, TN], 1:num1:TN};
for j=1:3
    if bu(j)
        bcu(bn{j})=false;
    end
    if bv(j)
        bcv(bn{j})=false;
    end
end
bc=[bcu, bcv];
Ke=Ke(bc,bc);
Me=Me(bc,bc);

% Solve Eigenvalues (SE): 1 - all, 2,3- first a few
Ns=10; SE=3;
[d, V, NT]=SolveEig(Ke, Me, Ns, SE);
[sd, I]=sort(d); 
Gs=E/(2*(1+nu));
sd=real(sqrt(sd'*rho/Gs));

% Plot mode shapes
k=4; dd=0.1;
Vt=zeros(2*TN, 1);
Vt(bc)=V(:, I(k));
Uv=Vt(1:TN); Vv=Vt(TN+(1:TN));
Uv=G*Uv; Vv=G*Vv; 
X=squeeze(pnts(1,:,:));
Y=squeeze(pnts(2,:,:));
Uv=reshape(Uv, m, n); 
Vv=reshape(Vv, m, n); 
figure; surf(X,Y,Uv); title('U');
figure; surf(X,Y,Vv); title('V');

mv=max(max(abs([Uv,Vv])));
Uv=Uv/mv; 
Vv=Vv/mv; 
Xd=X+dd*Uv;
Yd=Y+dd*Vv;
figure; hold on;
plot(Xd, Yd);
plot(Xd', Yd');
plot(X(:,[1,end]), Y(:,[1,end]), 'r-.');
plot(X([1,end],:)', Y([1,end],:)', 'r-.');
axis equal;
axis off;







