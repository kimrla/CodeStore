function crv=spiralcylind(a, b, sang, eang, n)

% Create a interpolated cylindrical spiral curve
% 
% Calling Sequences:
% 
%     crv=spiralcylind(a, b, sang, eang)
%     crv=spiralcylind(a, b, sang, eang, n)
%
%  Input: 
% 
%    a - the radius of the cylinder
%    b - the lead of the spiral curve
%    sang - the start angle
%    eang - the end angle
%    n - the number of points in each quarter of a circle
%         used to interplate the spiral curve
%         default value if 6 points in each quarter of a circle
%

if nargin==4
    n=12*round(abs(eang-sang)/(pi/2));
else
    n=n*round(abs(eang-sang)/(pi/2));
end

t=linspace(sang, eang, n);
x=a*cos(t);
y=a*sin(t);
z=b*t;
w=ones(size(x));
knots=linspace(0, 1, n);
crv=nrlmake([x; y; z; w], knots, [0, 1], 3);

%% Demo
% % Constants of the cylindrical spiral curve
% a=1; b=0.1; t1=0; t2=9*pi;
% 
% % Get the curve
% crv=spiralcylind(a, b, t1, t2);
% crvs(1)=crv;
% nrlplot(crv, 300);
% title('A cylindrical spiral curve');



