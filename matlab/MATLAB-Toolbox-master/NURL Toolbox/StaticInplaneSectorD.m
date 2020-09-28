clear; clc;

% Order of the NURLS basis (order)
% The number of nodes (nodes)
% Young's modulus (E), Poisson's ratio (nu)
% Density (rho) and forces (p1, p2)
E=72e9; nu=0.3; rho=2700; 
E=E/(1-nu^2); nu=nu/(1-nu); % Plane strain
order=5*[1, 1]; nodes=10*[2, 1];
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

% Stiffness matrix
C=E/(1-nu^2); nu1=(1-nu)/2;
G11=Gx'*CJ*Gx;
G22=Gy'*CJ*Gy;
G12=Gx'*CJ*Gy;
Ke=C*[G11+nu1*G22, nu*G12+nu1*G12';
           nu*G12'+nu1*G12, G22+nu1*G11];

% Force vector
m=length(s); n=length(t); 
[pnt, dp]=nrldeval(srf, {s, t}); 
X=squeeze(pnt(1,:,:)); Y=squeeze(pnt(2,:,:)); 
Xs=squeeze(dp{1}(1,:,:)); Ys=squeeze(dp{1}(2,:,:)); 
Xt=squeeze(dp{2}(1,:,:)); Yt=squeeze(dp{2}(2,:,:)); 
G1 = nrlgintvdeval (srf.order(1), srf.coefs(4,:,1), srf.knots{1}, s);
x1=X(:,1); y1=Y(:,1);
s1=Xs(:,1); s2=Ys(:,1); s3=0*s1;
ss1=[s1, s2, s3]; 
dn1=zeros(size(ss1)); dn1(:, 3)=1; 
G2 = nrlgintvdeval (srf.order(1), srf.coefs(4,:,end), srf.knots{1}, s);
x2=X(:,end); y2=Y(:,end);
s1=Xs(:,end); s2=Ys(:,end); s3=0*s1;
ss2=[s1, s2, s3]; 
dn2=zeros(size(ss2)); dn2(:, 3)=-1; 
dt1=cross(dn1, ss1, 2); 
dt2=cross(dn2, ss2, 2); 
Pu1=p1*G1*diag(Cs0)*dt1(:,1);
Pv1=p1*G1*diag(Cs0)*dt1(:,2);
Pu2=p2*G2*diag(Cs0)*dt2(:,1);
Pv2=p2*G2*diag(Cs0)*dt2(:,2);
Pn1=sqrt(Pu1.^2+Pv1.^2);
Pn2=sqrt(Pu2.^2+Pv2.^2);
B1=bcnodes([0,1,0,0], nodes(1),  nodes(2)); 
B2=bcnodes([0,0,0,1], nodes(1),  nodes(2)); 
Pu=0*G11(:,1); Pv=Pu;
Pu(B1)=Pu1(:); Pv(B1)=Pv1(:);
Pu(B2)=Pu2(:); Pv(B2)=Pv2(:);
Fe=[Pu; Pv];

% Apply boundary conditions
Bu=bcnodes(bu, nodes(1),  nodes(2)); 
Bv=bcnodes(bv, nodes(1),  nodes(2)); 
bc=~[Bu, Bv]; 
Ke=Ke(bc, bc); 
Fe=Fe(bc);

% Solve static deformation
Us=Ke\Fe;

% Draw static modes
Gs=E/(2*(1+nu));
NT=prod(nodes);
Ud=zeros(2*NT, 1);
Ud(bc)=Us(:);
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







