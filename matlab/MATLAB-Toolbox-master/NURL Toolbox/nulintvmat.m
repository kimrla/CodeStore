function G=nulintvmat(t, knots, order, der)

% Weighting coefficient matrix of Galerkin interpolation or its derivatives at an interval 
%   input : 
%               t - the evalued point 
%               knots - knots vector
%               order - the order of the Langrange basis
%               der - the order of derivatives
%  output :   
%  
%    G : Galerkin interpolation matrix or its derivatives
%

% Check order
if order>7
    order=7;
end

if isempty(t)
    G=[];
    return;
end

% Get span index
[index, iu]=spanweak(knots, t, order);

% Get the points
G=zeros(length(knots), length(t)); 
[n, ~]=size(iu); 
for i=1:n
    pp=index(i)+1:index(i+1);
    G(iu(i,:), pp)=nurlbasis(t(pp), knots(iu(i,:)), order, der); 
end

%% ! Test - basis
% % Oder of basis (order), derivatives (der), 
% % the basis number (k), and the span index (i)
% order=3; der=0; k=6; i=6;
% 
% % knots (N), number of sampling points (m)
% n=10; m=1001;
% 
% % Get vectors
% u=linspace(0, 1, n);
% t=linspace(0, 1, m);
% 
% % Get derivative by weighting matrix
% G=nulintvmat(t, u, order, der);
% 
% % Get span index
% [index, iu]=spanweak(u, t, order);
% 
% % Plot a few basis
% figure; plot(t, G([1,5,9],:));
% hold on;
% plot(u, 0*u, 'ro');
% q=length(index)-1;
% for i=1:q
%     pp=index(i)+1:index(i+1);
%     plot(t(pp), 0*t(pp), 'LineWidth', 2);
% end
% plot(u(3:end-1), 0*u(3:end-1), 'r*');
% title('A few NURL basis and the searching method');
% 
% % Plot all basis
% figure; plot(t, G);
% hold on;
% plot(u, 0*u, 'ro');
% title('All NURL basis');

%% Test - interpolation
% % A sine function and its first derivatives
% fun = @(t) sin(t);
% dfun =@(t) cos(t);
% 
% % Oder of basis (order) and derivatives (der)
% order=5; der=1;
% 
% % Knots (N), number of sampling points (m)
% n=10; m=10;
% 
% % Get knots (u) and interpolation nodes (t)
% u=linspace(-1, 1, n);
% t=linspace(-1, 1, m);
% 
% % Exact values
% y0=fun(u)';
% y=fun(t)'; dy=dfun(t)';
% 
% % Get derivative by weighting matrix
% G=nulintvmat(t, u, order, 0);
% A=nulintvmat(t, u, order, 1);
% fy=G'*y0;
% dfy=A'*y0;
% 
% % Plot values and errors
% figure; hold on;
% plot(t, y); plot(t, fy,'ro');
% legend('Exact', 'NURL', 'Location', 'northwest');
% title('Function values');
% 
% figure; hold on;
% plot(t, dy); plot(t, dfy,'ro'); 
% legend('Exact', 'NURL');
% title('First derivatives');
% 
% figure; plot(t, y-fy); title('y-fy');
% title('Error of function values');
% 
% figure; plot(t, dy-dfy);  title('dy-dfy');
% title('Error of  first derivatives');




