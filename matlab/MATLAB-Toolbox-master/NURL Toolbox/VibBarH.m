clear; clc;

% Oder of basis (order) and derivatives (der)
order=3; der=1;

% Length of the bar (L) and the number 
%   of sampling points (n)
L=1; n=200; 

% Get knots vectors and intergration nodes
u=linspace(0, L, n);
[x, C]=GaussLobattoQ(u, 2*order+2);
G=nuhintvmat(x, u, order, 0); 
A=nuhintvmat(x, u, order, 1); 

% Stiffness and mass matrices
C=repmat(C, 2*n, 1);
Ke=(A.*C)*A';
Me=(G.*C)*G';

% Applying boundary conditions
q=[2:2*n-2, 2*n];
Ke=Ke(q,q);
Me=Me(q,q);

% Solve Eigenvalues (SE): 1 - all, 2- first a few
Ns=50; SE=1;
[d, V, TN]=SolveEig(Ke, Me, Ns, SE);
[sd, I]=sort(d); 
sd=real(sd.^(1/2))/pi; 

% Draw mode shape
k=2;
Vt=[0; V(1:end-1,I(k)); 0; V(end, I(k))];
figure; plot(x, G'*Vt);

nn=length(sd);
figure; plot((1:nn)/nn, (sd-(1:nn)')./(1:nn)');





