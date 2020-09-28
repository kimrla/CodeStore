function dy=mfdm(x, y)

% Get derivatives of a function using midpoint finite difference method
% 
% Input: 
% 
%     x, y - the coordinates of the curve
% 
%  Output: 
% 
%    dy - the derivatives (dy/dx)
%

n=length(x); dy=zeros(size(y)); 
for i=2:n-1
    knots=[x(i-1), x(i), x(i+1)];
    dL=nurls(x(i), knots);
    dy(i)=dL(1)*y(i-1)+dL(2)*y(i)+dL(3)*y(i+1); 
end

i=1; knots=[x(i), x(i+1), x(i+2)];
dL=nurls(x(i), knots);
dy(1)=dL(1)*y(1)+dL(2)*y(2)+dL(3)*y(3); 

i=n; knots=[x(i-2), x(i-1), x(i)];
dL=nurls(x(i), knots);
dy(n)=dL(1)*y(n-2)+dL(2)*y(n-1)+dL(3)*y(n); 

function dL=nurls(x, knots)

x1=knots(1); x2=knots(2); x3=knots(3); 
dL=[-(x2 - 2.*x + x3)/((x1 - x2).*(x1 - x3))
        (x1 - 2.*x + x3)/((x1 - x2).*(x2 - x3))
       -(x1 - 2.*x + x2)/((x1 - x3).*(x2 - x3)) ];


