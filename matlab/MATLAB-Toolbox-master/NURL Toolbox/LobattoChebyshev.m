function x=LobattoChebyshev(a, b, N)

% Generate Lobatto Chebyshev nodes
% 
%  Input:
%    a, b - the span of the interval
%    N - the number of nodes required
%
%  Output:
%    x - Lobatto Chebyshev nodes
%

x=zeros(N,1); 
for j=1:N
    nd=(j-1)/(N-1); h=0.5*(1-cos(nd*pi));
    x(j)=a+h*(b-a);
end