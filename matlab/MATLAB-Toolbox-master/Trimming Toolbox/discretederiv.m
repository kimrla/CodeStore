function [dF,dF2]=discretederiv(f,x1,C)
% Use difference quotient to replace to obtain the derivative, only
% suitable to 1-D or 2-D.

% Input:
%   f: Original function handle.
%   x1: Par-coords of the point which needs to calculate the derivative at
%       this point.
%   c: Constant to confirm the step length, default is 0.01;
% Output:
%   dF: 1st derivative (for equation) or the Jacobi mattix (for equations).
%   dF2: 2nd derivative (for equation) or the Hessan mattix (for equations).

% dF=[Su,Sv], size(Su)=[3,1], size(dF)=[3,2], size(dF2)=[3,4]=[Suu,Suv,Suv,Svv]   
c=0.01;
if nargin==3
    c=C;
end
n=numel(x1);
F=f(x1); % size(F)=[3,nump], size(x0)=[1,1] or [1,2]
h=c.*norm(F).*ones(size(x1));% Step, size(h)=size(x), h~=0   
x2=x1+h;
x0=x1-h;
pp=find(x2>1);
h(pp)=-h(pp);
x2=x1+h;% Use backward difference quotient, if exceeds the domain (>1), 
% then use forward difference quotient.
H=repmat(h,[3,1]); % size(H)=[3,2] or [3,1]
if n==2 % Equations 
    dF=([f([x2(1),x1(2)]),f([x1(1),x2(2)])]-[F,F])./H; % Jacobi matrix
    if nargout==2 % Hessen matrix: [Suu,Suv,Svu,Svv]
        df=([F,F]-[f([x0(1),x1(2)]),f([x1(1),x0(2)])])./H;
        dF2(:,1)=(dF(:,1)-df(:,1))./H(:,1);
        dF2(:,4)=(dF(:,2)-df(:,2))./H(:,2);              
        tem1=((f(x2)-f([x1(1),x2(2)]))./H(:,1)-dF(:,1))./H(:,2);
        tem2=((f([x2(1),x1(2)])-f([x2(1),x0(2)]))./H(:,2)-dF(:,2))./H(:,1);
        dF2(:,2)=(tem1+tem2)/2;
        dF2(:,3)=dF2(:,2); % NEED TO CHECK THE FORMULA
    end
elseif n==1 % Equation
    dF=(f(x2)-f(x1))/h;
    dF2=(f(x2)+f(x0)-2*f(x1))/(2*h*h);
end
        





end