function dy=tfdm(x, y)

% Get derivatives of a function using two-point finite difference method
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
for i=1:n-1
    knots=[x(i), x(i+1)];
    dL=nurls(x(i), knots);
    dy(i)=dL(1)*y(i)+dL(2)*y(i+1); 
end

dy(n)=dy(n-1); 

function dL=nurls(x, knots)

x1=knots(1); x2=knots(2); 
dL1=[  1/(x1 - x2)
          -1/(x1 - x2)];
n=2; dL=ones(n, length(x));
for i=1:n
    dL(i,:)=dL1(i).*dL(i,:);
end


