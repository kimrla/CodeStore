function G=nuhintvmat(t, knots, order, der)

% Weighting coefficient matrix of Galerkin interpolation or its derivatives at an interval for NUH basis 
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

t=checktt(t);
knots=checktt(knots);

% Get span index
[index, iu]=spanweak(knots, t, order);

% Get the weights
G=zeros(2*length(knots), length(t)); 
[n, ~]=size(iu); 
for i=1:n
    pp=index(i)+1:index(i+1);
    qq=2*iu(i,1)-1:2*iu(i,end);
    G(qq, pp)=nurhbasis(t(pp), knots(iu(i,:)), order, der); 
end


%% Test - basis
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
% G=nuhintvmat(t, u, order, der);
% 
% % Get span index
% [index, iu]=spanweak(u, t, order);
% 
% % Plot a few basis
% figure; plot(t, G([2:2:20],:));
% hold on;
% plot(u, 0*u, 'ro');
% q=length(index)-1;
% for i=1:q
%     pp=index(i)+1:index(i+1);
%     plot(t(pp), 0*t(pp), 'LineWidth', 2);
% end
% plot(u(3:end-2), 0*u(3:end-2), 'r*');
% title('A few NURH basis and the searching method');
% 
% % Plot all basis
% figure; plot(t, G);
% hold on;
% plot(u, 0*u, 'ro');
% title('All NURH basis');

%% Demo - Interpolation
% % A sine function and its derivatives
% fun = @(t) sin(t);
% dfun =@(t) cos(t);
% ddfun =@(t) -sin(t);
% dddfun =@(t) -cos(t);
% 
% % Oder of basis (order) and derivatives (der)
% order=2; der=0;
% 
% % knots (N), number of sampling points (m)
% n=10; m=21;
% 
% % Get knots (u) and nodes (t)
% u=linspace(-1, 1, n);
% t=linspace(-1, 1, m);
% 
% % Exact values
% y0=fun(u); dy0=dfun(u); 
% y=fun(t); dy=dfun(t); ddy=ddfun(t); d3y=dddfun(t); 
% 
% % Get derivative by weighting matrix
% G=nuhintvmat(t, u, order, der);
% D1=nuhintvmat(t, u, order, der+1);
% D2=nuhintvmat(t, u, order, der+2);
% D3=nuhintvmat(t, u, order, der+3);
% 
% n=length(u); 
% p=1:2:2*n-1; q=2:2:2*n; 
% yr=zeros(2*n,1);
% yr(p)=y0; yr(q)=dy0; 
% yv=G'*yr;
% ydv=D1'*yr;
% yd2v=D2'*yr;
% yd3v=D3'*yr;
% 
% % Plot results
% figure; hold on;
% plot(t, y); plot(t, yv, 'ro'); 
% title('Function values');
% 
% figure; plot(t, y-yv'); 
% title('Error of function values');
% 
% figure; hold on;
% plot(t, dy); plot(t, ydv, 'ro'); 
% title('First derivatives');
% 
% figure; plot(t, dy-ydv'); 
% title('Error of first derivatives');
% 
% figure; hold on;
% plot(t, ddy); plot(t, yd2v, 'ro'); 
% title('Second derivatives');
% 
% figure; plot(t, ddy-yd2v'); 
% title('Error of second derivatives');
% 
% figure; hold on;
% plot(t, d3y); plot(t, yd3v, 'ro'); 
% title('Third derivatives');
% 
% figure; plot(t, d3y-yd3v'); 
% title('Error of third derivatives');






