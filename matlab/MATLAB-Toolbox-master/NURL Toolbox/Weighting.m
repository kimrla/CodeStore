function A = Weighting(x, i)

% Get weighting coefficient matrix of DQM
% 
% Calling Sequences:
% 
%     A = WeightingI(x)
% 
%     Ai = WeightingI(x, i)
%   
%  Input:
% 
%    x - knots vector
% 
%    i - index of a node
%
%  Output:
% 
%     A - weighting coefficient matrix
%  
%    Ai - weighting coefficients of node i
% 

N=length(x); 
if nargin==2
    A=zeros(1, N);
    for j=1:N
        if j==i
            pp=true(N,1); pp(i)=false; 
            A(j)=sum(1./(x(i)-x(pp)));
        else
            pp=true(N,1); qq=pp;
            pp([i,j])=false;
            qq(j)=false;
            A(j)=prod(x(i)-x(pp))/prod(x(j)-x(qq));
        end
    end
elseif nargin==1
    A=zeros(N);
    for i=1:N
        A(i, :)=Weighting(x, i);
    end
end





