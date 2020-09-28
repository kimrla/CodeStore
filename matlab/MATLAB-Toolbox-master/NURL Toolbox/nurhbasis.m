function dnH=nurhbasis(x, knots, p, der)

% Get the basis of non-uniform Hermite basis
% 
%   input :  x - the evalued point
%               knots - knot vector of each basis
%               order - the order of the Langrange basis
%               der - the order of derivatives
% 
%  output : dnL - a matrix includes the nurhs basis or 
%                         their derivatives
%

dnB1=zeros(p+1, length(x)); 
dnB2=dnB1;
dnL=cell(der+1, 1);
for i=1:der+1
    dnL{i}=nurlbasis(x, knots, p, i-1); 
end
if der==0
    dnB2=dnL{1}.*dnL{1};
elseif der==1
    dnB1=dnL{1}.*dnL{1};
    dnB2=2*dnL{1}.*dnL{2};
else
    for i=1 : der
        if i<der
            dnB1=dnB1+2*nchoosek(der-2, i-1)*dnL{i}.*dnL{der-i+1};
        end
        dnB2=dnB2+2*nchoosek(der-1, i-1)*dnL{i}.*dnL{der-i+2};
    end
end
dnH=zeros(2*(p+1), length(x));
dnLi=nurlbasis(knots, knots, p, 1); 
for i=1:p+1
    a=-2*dnLi(i,i); b=2*knots(i)*dnLi(i,i)+1; 
    dnH(2*i-1,:)=der*a*dnB1(i,:)+(a*x+b).*dnB2(i,:); 
    a=1; b=-knots(i); 
    dnH(2*i,:)=der*a*dnB1(i,:)+(a*x+b).*dnB2(i,:); 
end


%% Demo
% p=3; % Order of basis
% der=0; % Order of derivatives
% n=100; % Number of nodes
% a=-1; b=1; % Span of the function
% 
% k=4; % Basis index number
% x=LobattoChebyshev(a, b, n); x=x';
% knots=linspace(a, b, p+1);
% 
% dnH=nurhbasis(x, knots, p, der); 
% ddnH=nurhbasis(x, knots, p, der+1);
% 
% A=Weighting(x);
% ddH=dnH*A';
% 
% figure; hold on;
% plot(x,dnH(k,:)); 
% plot(knots,0*knots,'ro');
% title('The basis');
% 
% figure; hold on;
% plot(x, ddnH(k,:)); 
% plot(x, ddH(k,:), 'k*'); 
% plot(knots, 0*knots, 'ro');
% title('First derivatives');






