clear; clc;

% Order of the NURLS basis (order)
% The number of nodes (nodes)
% Young's modulus (E), Poisson's ratio (nu)
% Density (rho) and forces (p1, p2)
E=72e9; nu=0.3; rho=2700; 
E=E/(1-nu^2); nu=nu/(1-nu); % Plane strain
order=5*[1, 1]; nodes=15*[2, 1];
p1=-100; p2=-10;

% Boundary conditions (bu, bv)
% 1 - clamped, 0 - free
bu=[0 0 1 0]; bv=[1 0 0 0]; 

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
[s, Cs0]=GaussLobattoQ(srf.knots{1}, srf.order(1)+2);
[t, Ct0]=GaussLobattoQ(srf.knots{2}, srf.order(2)+2);
[Ct, Cs]=meshgrid(Ct0, Cs0);
[intpi, jac, Jacb] = nrlsrfgintvdeval (srf.order, srf.coefs, srf.knots, {s, t});
G=intpi'; Gx=jac{1}'; Gy=jac{2}'; 
CJ=Jacb.*Cs(:)'.*Ct(:)';
CJ=diag(CJ);

% Stiffness and mass matrices
C=E/(1-nu^2); nu1=(1-nu)/2;
G11=Gx'*CJ*Gx;
G22=Gy'*CJ*Gy;
G12=Gx'*CJ*Gy;
Ke=C*[G11+nu1*G22, nu*G12+nu1*G12';
           nu*G12'+nu1*G12, G22+nu1*G11];
Gm=G'*CJ*G; 
Zm=zeros(size(Gm)); 
Me=rho*[Gm, Zm; Zm, Gm]; 

% Static deformation
x=squeeze(srf.coefs(1,:,:)); 
y=squeeze(srf.coefs(2,:,:)); 
[Tx,Ty,Txy,Trr,Ttt,Trt,Ui,Vi,Ur]=ThickCylinderAnalytical(x,y,p1,p2);
Bs=bcnodes([0,1,0,1], nodes(1),  nodes(2)); 
uv0=[Ui(:); Vi(:)]; Bs0=[Bs, Bs];

% Apply boundary conditions
Bu=bcnodes(bu, nodes(1),  nodes(2)); 
Bv=bcnodes(bv, nodes(1),  nodes(2)); 
bc=~[Bu, Bv]; 
Ke=Ke(bc, bc); 
Me=Me(bc, bc); 
Bs=Bs0(bc);
uv=uv0(bc);

% Solve Eigenvalues (SE): 1 - all, 2, 3- first a few
Ns=10; SE=3;
[d, V, TN]=SolveEig(Ke, Me, Ns, SE);
[sd, I]=sort(d); 
Gs=E/(2*(1+nu));
sd=real(sd.^(1/2))*sqrt(rho/Gs); 

% Solve static deformation
Bk=~Bs;
Ks=Ke(Bk, Bk);
Fs=-Ke(Bk, Bs)*uv(Bs);
Us=Ks\Fs;

% Draw vibration mode shapes
k=1;
m=length(s); n=length(t); 
[pnt, dp]=nrldeval(srf, {s, t});
X=squeeze(pnt(1,:,:)); 
Y=squeeze(pnt(2,:,:)); 
NT=prod(nodes);
Vt=zeros(2*NT, 1);
Vt(bc)=V(:, I(k));
Ut=Vt(1:NT);
Vt=Vt(NT+(1:NT));
Ut=reshape(G*Ut, m, n); 
Vt=reshape(G*Vt, m, n); 
figure; surf(X, Y, Ut); title('U'); 
figure; surf(X, Y, Vt); title('V'); 

% Draw static deformation
Uk=zeros(size(Bk))'; 
Uk(Bk)=Us(:); 
Uk(Bs)=uv(Bs);
Ud=zeros(2*NT, 1); 
Ud(bc)=Uk(:);
Ut=Ud(1:NT);
Vt=Ud(NT+(1:NT));
Tx=C*(Gx*Ut+nu*Gy*Vt);
Ty=C*(Gy*Vt+nu*Gx*Ut);
Txy=Gs*(Gx*Vt+Gy*Ut);
Ux=reshape(G*Ut, m, n); 
Uy=reshape(G*Vt, m, n); 
Tx=reshape(Tx, m, n); 
Ty=reshape(Ty, m, n); 
Txy=reshape(Txy, m, n); 
figure; surf(X, Y, Ux); title('U'); 
figure; surf(X, Y, Uy); title('V'); 
figure; surf(X, Y, Tx); title('Tx'); 
figure; surf(X, Y, Ty); title('Ty'); 
figure; surf(X, Y, Txy); title('Txy'); 

R=sqrt(X.^2+Y.^2);
cost = X./R; sint=Y./R;
Ur=Ux.*cost+Uy.*sint;
Ut=Uy.*cost-Ux.*sint;
Tr=Tx.*cost.^2+Ty.*sint.^2+2*Txy.*sint.*cost;
Tt=Tx.*sint.^2+Ty.*cost.^2-2*Txy.*sint.*cost;
Trt=Txy.*(cost.^2-sint.^2)+(Ty-Tx).*sint.*cost;

figure; surf(X, Y, Ur); title('Ur'); 
figure; surf(X, Y, Ut); title('Ut'); 
figure; surf(X, Y, Tr); title('Tr'); 
figure; surf(X, Y, Tt); title('Tt'); 
figure; surf(X, Y, Trt); title('Trt'); 


