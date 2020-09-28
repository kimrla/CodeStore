clear; clc;

% Oder of basis (order) and derivatives (der)
order=5; der=1;

% Length of the bar (L) and the number 
%   of sampling points (n)
L=1; n=300; 

% Get knots vectors and intergration nodes
u=linspace(0, L, n);
span1=u(1 : end-1); 
span2=u(2 : end); 
s=sort([u, (span1+span2)/2]); 
[x, C]=GaussLobattoQ(u, order+2);
G=nulintvmat(x, u, order, 0);
A=nulintvmat(x, u, order, 1);

% Stiffness and mass matrices
C=repmat(C, n, 1);
Ke=(A.*C)*A';
Me=(G.*C)*G';

% Applying boundary conditions
q=2:n-1;
Ke=Ke(q,q);
Me=Me(q,q);

% Solve Eigenvalues (SE): 1 - all, 2- first a few
Ns=50; SE=1;
[d, V, TN]=SolveEig(Ke, Me, Ns, SE);
[sd, I]=sort(d); 
sd=real(sd.^(1/2))/pi; 

% Draw mode shape
k=2;
Vt=[0; V(:,I(k)); 0];
figure; plot(x, G'*Vt);

nn=length(sd);
figure; plot((1:nn)/nn, (sd-(1:nn)')./(1:nn)');





