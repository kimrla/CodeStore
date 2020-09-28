function [d, V, N]=SolveEig(K, M, N, SE)

[TN,~]=size(K);
switch SE
    case 1
        N=TN;
        [V, Di]=eig(K, M);
        d=zeros(N,1);
        for j=1:N
            d(j,1)=Di(j,j);
        end
    case 2
        % Cholesky factorization
        Rc = chol(M, 'upper'); Ri=inv(Rc); 
        K=Ri'*K*Ri;       
        
        % Solve Eigenvalues
        [V,Di]=eigs(K, N, 'sm'); d=zeros(N,1);
        for j=1:N
            d(j,1)=Di(j,j);
        end
        V=Ri*V;
    case 3
        [V,Di]=eigs(K, M, N, 'sm'); d=zeros(N,1);
        for j=1:N
            d(j,1)=Di(j,j);
        end
end






