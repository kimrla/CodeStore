% Vibration analysis of sectorial thin plate by NURL
clear; clc;

% Order of the NURLS basis (order)
% The number of nodes (nodes)
% Young's modulus (E), Poisson's ratio (nu), density (rho) and thickness (h)
E=200e9; nu=0.3; rho=2700; h=0.05;
order=5*[1, 1]; nodes=12*[2, 1];

% Boundary conditions (bw)
% 1 - clamped, 0 - free
bc=[1 1 1 1]; 

% Create a sector
crv1 = nrlcirc(0.5, [0, 0], 0, pi/4);
crv2 = nrlcirc(1, [0, 0], 0, pi/4);
srf=nrlruled(crv1, crv2);
srf = nrldegelev(srf, order-srf.order);
iknts=nodes-srf.number;
srf=nrlkntins(srf, {iknts(1), iknts(2)});
figure; nrlplot(srf, [100, 100], 'ctrl');
view(2); axis equal;

% Get first derivatives and weights of derivatives in Cartesian coordinates
[s, Cs]=GaussLobattoQ(srf.knots{1}, srf.order(1)+5);
[t, Ct]=GaussLobattoQ(srf.knots{2}, srf.order(2)+5);
[Ct, Cs]=meshgrid(Ct, Cs);
[intpi, ~, geom, hess] = nrlsrfgintvdeval (srf.order, srf.coefs, srf.knots, {s, t});
G=intpi'; Gxx=hess{1,1}'; Gyy=hess{2,2}'; Gxy=hess{1,2}'; 
CJ=geom{1}.*Cs(:)'.*Ct(:)';
CJ=repmat(CJ', 1, size(intpi, 1));

% Stiffness and mass matrices
D=E*h^3/(12*(1-nu^2)); 
Ke=D*(Gxx'*(CJ.*Gxx)+Gyy'*(CJ.*Gyy)+nu*(Gxx'*(CJ.*Gyy)+Gyy'*(CJ.*Gxx))+2*(1-nu)*Gxy'*(CJ.*Gxy));
Me=rho*h*G'*(CJ.*G);

% Apply boundary conditions
bc=~bcnodes(bc, nodes(1),  nodes(2)); 
Ke=Ke(bc, bc); 
Me=Me(bc, bc); 

% Solve Eigenvalues (SE): 1 - all, 2, 3- first a few
Ns=30; SE=3;
[d, V, TN]=SolveEig(Ke, Me, Ns, SE);
[sd, I]=sort(d); 
sd=real(sd.^(1/2))*sqrt(rho*h/D); 

% Draw mode shapes
k=2;
pnt=nrleval(srf, {s, t});
m=length(s); n=length(t); 
NT=prod(nodes);
W=zeros(NT, 1);
W(bc)=V(:, I(k));
W=reshape(G*W, m, n); 
X=squeeze(pnt(1,:,:));
Y=squeeze(pnt(2,:,:));
figure; surf(X, Y, W); 






