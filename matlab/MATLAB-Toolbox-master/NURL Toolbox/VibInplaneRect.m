clear; clc;

% Order of the NURL basis (order)
% The number of nodes (nodes)
% Young's modulus (E), Poisson's ratio (nu) and density (rho)
E=200e9; nu=0.3; rho=2700; 
order=5*[1, 1]; nodes=20*[1, 1];

% Boundary conditions (bu, bv)
% 1 - clamped, 0 - free
bu=[0 1 0 1]; bv=[1 0 1 0]; 

% Create a sector
crv1 = nrlline([0,0], [1,0]);
crv2 = nrlline([0,1], [1,1]);
srf=nrlruled(crv1, crv2);
srf = nrldegelev(srf, order-srf.order);
iknts=nodes-srf.number;
srf=nrlkntins(srf, {iknts(1), iknts(2)});
nodes=srf.number;
figure; nrlplot(srf, [100, 100], 'ctrl');
view(2); axis equal;

% Get first derivatives and weights of derivatives in Cartesian coordinates
[s, Cs]=GaussLobattoQ(srf.knots{1}, srf.order(1)+2);
[t, Ct]=GaussLobattoQ(srf.knots{2}, srf.order(2)+2);
[Ct, Cs]=meshgrid(Ct, Cs);
[intpi, jac, geom] = nrlsrfgintvdeval (srf.order, srf.coefs, srf.knots, {s, t});
G=intpi'; Gx=jac{1}'; Gy=jac{2}'; 
dpxi=jac{1}; dpyi=jac{2}; 
CJ=geom{1}.*Cs(:)'.*Ct(:)';
CJ=repmat(CJ', 1, size(intpi, 1));

% Stiffness and mass matrices
C=E/(1-nu^2); nu1=(1-nu)/2;
G11=Gx'*(CJ.*Gx);
G22=Gy'*(CJ.*Gy);
G12=Gx'*(CJ.*Gy);
Ke=C*[G11+nu1*G22, nu*G12+nu1*G12';
           nu*G12'+nu1*G12, G22+nu1*G11];
Gm=G'*(CJ.*G); 
Zm=zeros(size(Gm)); 
Me=rho*[Gm, Zm; Zm, Gm]; 

% Apply boundary conditions
Bu=bcnodes(bu, nodes(1),  nodes(2)); 
Bv=bcnodes(bv, nodes(1),  nodes(2)); 
bc=~[Bu, Bv]; 
Ke=Ke(bc, bc); 
Me=Me(bc, bc); 

% Solve Eigenvalues (SE): 1 - all, 2, 3- first a few
Ns=30; SE=1;
[d, V, TN]=SolveEig(Ke, Me, Ns, SE);
[sd, I]=sort(d); 
Gs=E/(2*(1+nu));
sd=real(sd.^(1/2))*sqrt(rho/Gs)/pi; 

% Draw mode shapes
k=1;
pnt=nrleval(srf, {s, t});
m=length(s); n=length(t); 
NT=prod(nodes);
Vt=zeros(2*NT, 1);
Vt(bc)=V(:, I(k));
Ut=Vt(1:NT);
Vt=Vt(NT+(1:NT));
Ut=reshape(G*Ut, m, n); 
Vt=reshape(G*Vt, m, n); 
X=squeeze(pnt(1,:,:));
Y=squeeze(pnt(2,:,:));
figure; surf(X, Y, Ut); title('U');
figure; surf(X, Y, Vt); title('V');

% Plot comparision curves
load exs1 esdd;
esdd=esdd(1:TN);
figure; plot((1:TN)/TN,(sd(:)-esdd)./esdd)

sdd=sd';





