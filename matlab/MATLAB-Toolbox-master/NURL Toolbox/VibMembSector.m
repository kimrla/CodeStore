% Vibration of a membrane sector
clear; clc;

% Order of the NURLS basis (order)
% The number of nodes (nodes)
% Surface tension (S) and surface density (rho)
order=5*[1, 1]; nodes=10*[2, 1];
S=1; rho=1; 

% Boundary conditions (bw): 1 - clamped, 0 - free
bw=[1 1 1 1]; 

% Create a sector
crv1 = nrlcirc(0.5, [0, 0], 0, pi/2);
crv2 = nrlcirc(1, [0, 0], 0, pi/2);
srf=nrlruled(crv1, crv2);
srf = nrldegelev(srf, order-srf.order);
iknts=nodes-srf.number;
srf=nrlkntins(srf, {iknts(1), iknts(2)});
figure; nrlplot(srf, [100, 100], 'ctrl');
view(2); axis equal;

% Get first derivatives and weights of derivatives in Cartesian coordinates
[s, Cs]=GaussLobattoQ(srf.knots{1}, srf.order(1)+2);
[t, Ct]=GaussLobattoQ(srf.knots{2}, srf.order(2)+2);
[Ct, Cs]=meshgrid(Ct, Cs);
[intpi, jac, geom] = nrlsrfgintvdeval (srf.order, srf.coefs, srf.knots, {s, t});
G=intpi'; Gx=jac{1}'; Gy=jac{2}'; 
CJ=geom{1}.*Cs(:)'.*Ct(:)';
CJ=repmat(CJ', 1, size(intpi, 1));

% Stiffness and mass matrices
Ke=S*(Gx'*(CJ.*Gx)+Gy'*(CJ.*Gy));
Me=rho*G'*(CJ.*G); 

% Apply boundary conditions
Bw=bcnodes(bw, nodes(1),  nodes(2)); 
bc=~Bw; 
Ke=Ke(bc, bc); 
Me=Me(bc, bc); 

% Solve Eigenvalues (SE): 1 - all, 2, 3- first a few
Ns=30; SE=3;
[d, V, TN]=SolveEig(Ke, Me, Ns, SE);
[sd, I]=sort(d); 
sd=real(sd.^(1/2)); 

% Draw mode shapes
k=6;
pnt=nrleval(srf, {s, t});
m=length(s); n=length(t); 
NT=prod(nodes);
W=zeros(NT, 1);
W(bc)=V(:, I(k));
W=reshape(G*W, m, n); 
X=squeeze(pnt(1,:,:));
Y=squeeze(pnt(2,:,:));
figure; hold on;
plot(X([1,end],:)', Y([1,end],:)', 'k');
plot(X(:,[1,end]), Y(:,[1,end]), 'k');
contour(X, Y, W, 20);
axis equal;


