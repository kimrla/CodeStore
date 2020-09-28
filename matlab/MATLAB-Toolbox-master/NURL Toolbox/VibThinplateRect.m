% Vibration analysis of rectangular thin plate by NUH
clear; clc; 

% Order of the NURLS basis (order)
% The number of nodes (nodes)
% Young's modulus (E), Poisson's ratio (nu), density (rho) and thickness (h)
E=200e9; nu=0.3; rho=2700; h=0.1; 
order=6*[1, 1]; nodes=8*[1, 1]; 

% Boundary conditions (bw, bwx, bwy)
% 1 - clamped, 0 - free
bw=[1 1 1 1]; bwx=[0 1 0 1]; 
bwy=[1 0 1 0]; bwxy=[0 0 0 0]; 

% Create a rectangular plate
crv1 = nrlline([0, 0], [1,0]);
crv2 = nrlline([0, 1], [1, 1]);
srf=nrlruled(crv1, crv2);
srf = nrldegelev(srf, order-srf.order);
iknts=nodes-srf.number;
srf=nrlkntins(srf, {iknts(1), iknts(2)});
figure; nrlplot(srf, [100, 100], 'ctrl');
view(2); axis equal;

% Get first derivatives and weights of derivatives in Cartesian coordinates
[s, Cs]=GaussLobattoQ(srf.knots{1}, 2*srf.order(1)+2);
[t, Ct]=GaussLobattoQ(srf.knots{2}, 2*srf.order(2)+2);
[Ct, Cs]=meshgrid(Ct, Cs);
[intpi, jac, geom, hess] = nuhsrfgintvdeval (srf, srf.order, srf.knots, {s, t});
G=intpi'; Gx=jac{1}'; Gy=jac{2}';
Gxx=hess{1,1}'; Gyy=hess{2,2}'; Gxy=hess{1,2}'; 
[pt, jc, hs] = nuhindex(nodes(1), nodes(2));
bn=[find(pt); find(jc{1}); find(jc{2}); find(hs)];
G=G(:, bn); Gx=Gx(:, bn); Gy=Gy(:, bn); 
Gxx=Gxx(:, bn); Gyy=Gyy(:, bn); Gxy=Gxy(:, bn); 
CJ=geom{1}.*Cs(:)'.*Ct(:)';
CJ=repmat(CJ', 1, length(bn));

% Stiffness and mass matrices
D=E*h^3/(12*(1-nu^2)); 
Ke=D*(Gxx'*(CJ.*Gxx)+Gyy'*(CJ.*Gyy)+nu*(Gxx'*(CJ.*Gyy)+Gyy'*(CJ.*Gxx))+2*(1-nu)*Gxy'*(CJ.*Gxy));
Me=rho*h*G'*(CJ.*G);

% Transform boundary DOFs
[bp, tbp, nbp]=boudvec(srf);
[Qx, Qn, bv, bt]=nuhsrftransmat(srf, tbp, nbp);
Ke=nuhsrftrans(Ke, bv, bt, Qx, Qn);
Me=nuhsrftrans(Me, bv, bt, Qx, Qn);

% Apply boundary conditions
bcw=~bcnodes(bw, nodes(1),  nodes(2)); 
bcwx=~bcnodes(bwx, nodes(1),  nodes(2)); 
bcwy=~bcnodes(bwy, nodes(1),  nodes(2)); 
bcwxy=~bcnodes(bwxy, nodes(1),  nodes(2)); 
bc=[bcw, bcwx, bcwy, bcwxy];
Ke=Ke(bc, bc); 
Me=Me(bc, bc); 

% Solve Eigenvalues (SE): 1 - all, 2, 3- first a few
Ns=30; SE=3; 
[d, V, TN]=SolveEig(Ke, Me, Ns, SE); 
[sd, I]=sort(d); 
sd=real(sd.^(1/2))*sqrt(rho*h/D)/pi^2; 
d1=1:3*nodes(1); d2=1:3*nodes(2); 
[d1, d2]=meshgrid(d1, d2); 
se=d1.^2+d2.^2; 
se=sort(se(:)); 

% Get mode shapes
k=4; 
m=length(s); n=length(t);  
NT=4*prod(nodes); 
W0=zeros(NT, 1); 
W0(bc)=V(:, I(k)); 
X=squeeze(geom{2}(1,:,:)); 
Y=squeeze(geom{2}(2,:,:)); 

% Transform back
W0=nuhsrftransbk(W0, bv, bt, Qx, Qn);

% Draw mode shapes
W=reshape(G*W0, m, n);  
figure; surf(X, Y, W); title('W'); 

Wx=reshape(Gx*W0, m, n); 
figure; surf(X, Y, Wx); title('Wx'); 

Wy=reshape(Gy*W0, m, n); 
figure; surf(X, Y, Wy); title('Wy'); 

figure; plot((1:TN)/TN, (sd-se(1:TN))./se(1:TN)); 





