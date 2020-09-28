function [index, iu]=spanweak(u, t, order)

% Get the spans and indexes of evaluation points
% 
%   input : u - the knots vector
%              t - evaluation points
%            order - order of the basis
% 
%   output :  iu - spans of evaluation points
%                indexes - indexes of the spans
%
% See also: nulintvmat

% Get span index
n=length(u); m=length(t);
iu=zeros(n-order, order+1);
iu(1:n-order, 1)=1:n-order;
for i=2:order+1
    iu(1:n-order, i)=iu(1:n-order, i-1)+1;
end
p=fix(order/2)+1;
index=zeros(n-order+1, 1);
for i=2:n-order
    pp=find(t>=u(i+p-1), 1)-1;
    if isempty(pp)
        index(i)=m;
    else
        index(i)=pp;
    end
end
index(end)=m;









