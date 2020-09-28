function crv=nrlanalytical(p, N, a1, a2, funx, funy, funz, funw)

% Create a curve using analytical functions
% 
% Calling Sequences:
% 
%     crv=nrlanalytical(p, N, a1, a2, funx, funy)
%     crv=nrlanalytical(p, N, a1, a2, funx, funy, funz)
%     crv=nrlanalytical(p, N, a1, a2, funx, funy, funz, funw)
%
%  Inputs :  
%                 p  :    the order of NURL basis
%                N  :    the number of control points
%               a1  :   the lower limit of the defined interval 
%               a2  :   the upper limit of the defined interval 
%            funx :  parametric function for x-coordinate
%            funy :  parametric function for y-coordinate 
%            funz :  parametric function for z-coordinate
%            funw :  parametric function of weights
% 
%  Example   :  a1=0; a2=2*pi; 
%             x0=0; y0=0; a=4;  b=3;
%             p=2; N=15;
%             funx = @(t) x0 + a*cos(t);
%             funy = @(t) y0 + b*sin(t);
%             crv=nrlanalytical(p, N, a1, a2, funx, funy);
% 

if nargin==6
    funz = @(t) 0*t;
    funw = @(t) ones(size(t));
elseif nargin==7
    funw = @(t) ones(size(t));
end

knots=linspace(0, 1, N); 
t=(a2-a1)*knots+a1; 
w=funw(t); 
x=funx(t)./w; 
y=funy(t)./w; 
z=funz(t)./w; 

crv=nrlmake([x; y; z; w], knots, [0, 1], p);


%% Demo - nul
% a1=0; a2=2*pi; 
% x0=0; y0=0; a=4;  b=3;
% p=2; N=15;
% funx = @(t) x0 + a*cos(t);
% funy = @(t) y0 + b*sin(t);
% crv=nrlanalytical(p, N, a1, a2, funx, funy);
% nrlplot(crv);

%% Demo - nurl
% p=2; N=3;
% a1=0; a2=1; 
% funx = @(t) (1-t.^2);
% funy = @(t) 2*t;
% funz = @(t) 0*t;
% funw = @(t) 1+t.^2;
% crv=nrlanalytical(p, N, a1, a2, funx, funy, funz, funw);
% figure; nrlctrlplot(crv);


