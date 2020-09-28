% Deriving the nurl basis by symbolic computation
clear; clc

n=5;

syms x 

% Get the variables and the coefficients
X=sym(zeros(n+1, 1));
xx=sym(zeros(1, n+1));
A=sym(zeros(n+1));
for i=1:n+1
    X(i)=sym(['x', num2str(i')]);
    xx(i)=x^(i-1);
    for j=1:n+1
        A(i, j)=X(i)^(j-1);
    end
end

% Get the Langrange basis
L=sym(zeros(n+1, 1));
for i=1:n+1
    b=zeros(n+1,1);
    b(i)=1;
    L(i)=simplify(xx*(A\b));
end

% The n-th order derivative of L
dnL=sym(zeros(n+1,n+1));
dnL(:,1)=L(:);
for i=2:n+1
    dnL(:,i)=simplify(diff(dnL(:,i-1), x));
end

% Vlaues at the interpolation points
vdnL=sym(zeros(n+1,n+1));
for i=1:n+1
    vdnL(:,i)=simplify(subs(dnL(:,i), x, X(i)));
end

dnL(:,1)



