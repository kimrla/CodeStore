function isl=isnrlline(crv, n)

% Check whether a nurl curve is a straight line
% 
%  Inputs: 
% 
%     crv - a nurls curve
% 
% Output: 
% 
%     isl - true (1) or false (0)
% 
%  Examples:
% 
%     lin=nrlline([0,5], [8,3]);
%     isl=isnrlline(lin);
% 

if nargin==1
    n=10;
end
pnt1=nrleval(crv, linspace(0, 1, n));   
pnt2=pnt1(:,2:n); pnt1(:,n)=[];
vect1=vecnorm(pnt2-pnt1);
vect2=vect1(:,2:n-1); vect1(:,n-1)=[];
C = cross(vect1, vect2);
nm=norm(C);

if nm<1e-10
    isl = 1;
else
    isl = 0;
end
        


