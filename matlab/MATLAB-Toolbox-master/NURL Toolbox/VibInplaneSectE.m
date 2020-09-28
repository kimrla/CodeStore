% Vibration analysis of sectorial plate by NURL as common FEM
clear; clc;

% Order of the NURLS basis (order)
% The number of nodes (nodes)
% Young's modulus (E), Poisson's ratio (nu) and density (rho)
E=200e9; nu=0.3; rho=2700; 
order=3*[1, 1]; ins=5*[2, 1];

% Boundary conditions (bu, bv)
% 1 - clamped, 0 - free
bu=[0 1 1 0]; bv=[1 1 0 0]; 

% Create a sector
crv1 = nrlcirc(0.5, [0, 0], 0, pi/2);
crv2 = nrlcirc(1, [0, 0], 0, pi/2);
srf=nrlruled(crv1, crv2);
srf = nrldegelev(srf, order-srf.order);
iu=linspace(0,1,ins(1)+1); 
iv=linspace(0,1,ins(2)+1); 
srf=nrlintins(srf, {iu(2:end-1), iv(2:end-1)});
figure; nrlplot(srf, [100, 100], 'ctrl');
view(2); axis equal;

% Node numbers of the element
num1=srf.number(1); 
num2=srf.number(2); 
TN=num1*num2;
nodes=zeros(num1, num2);
t=1;
for j=1:num2
    for i=1:num1
        nodes(i,j)=t;
        t=t+1;
    end
end

% Stiffness and mass matrcies of each element
em=length(srf.intervals{1})-1;
en=length(srf.intervals{2})-1;
order=srf.order;
knots={linspace(0,1,order(1)+1), linspace(0,1,order(2)+1)};
K=zeros(2*TN); M=K; 
elemat = cell(em*en, 1); 
ne=prod(order+1);
pp=zeros(1,2*ne);
pp(1:2:2*ne-1)=1:ne;
pp(2:2:2*ne)=ne+(1:ne);
for i=1:em
    p=(i-1)*order(1)+1:i*order(1)+1;
    for j=1:en
        q=(j-1)*order(2)+1:j*order(2)+1;
        coefsi=srf.coefs(:,p,q);
        Num=nodes(p,q);
        srfi=nrlmake(coefsi, knots, {[0,1],[0,1]}, order);
        [Ke, Me, pnt, G, Gx, Gy]=StiffMassInp(srfi, E, nu, rho);
        Ke=Ke(pp, pp); Me=Me(pp, pp);
        K=AssembleInP(K, Ke, Num(:));
        M=AssembleInP(M, Me, Num(:));
        k=(j-1)*em+i;
        elemat{k}.srf=srfi;
        elemat{k}.Num=Num;
        elemat{k}.pnt=pnt;
        elemat{k}.G=G;
        elemat{k}.Gx=Gx;
        elemat{k}.Gy=Gy;
    end
end

% Apply boundary conditions
ne=num1*num2;
pp=zeros(1,2*ne);
pp(1:2:2*ne-1)=1:ne;
pp(2:2:2*ne)=ne+(1:ne);
Bu=bcnodes(bu, num1,  num2); 
Bv=bcnodes(bv, num1,  num2); 
bc=~[Bu, Bv]; 
bc=bc(pp);
K=K(bc, bc); 
M=M(bc, bc); 

% Solve Eigenvalues (SE): 1 - all, 2, 3- first a few
Ns=10; SE=3;
[d, V, TN]=SolveEig(K, M, Ns, SE);
[sd, I]=sort(d); 
Gs=E/(2*(1+nu));
sd=real(sd'.^(1/2))*sqrt(rho/Gs); 





