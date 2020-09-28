function crv=nrl3pntsarc(pnt1, pnt2, pnt3)

% Create a nurl arc (<= pi) by three points in sequence
%  Inputs: 
% 
%      pnt1 - point 1
%      pnt2 - point 2
%      pnt3 - point 3
%
% Output:
% 
%      crv - a nurl circle
%

% Get the center of the triangle formed by the three points
center= planeline('center3ps', pnt1, pnt2, pnt3);

% Get the raduius
R=norm(pnt1-center);

% Get the vectors from center to points
v1=pnt1 - center; 
v2=pnt2 - center; 
v3=pnt3 - center; 

% Make the points be anticlockwise 
t1=[v1, 0]; t2=[v2, 0]; 
vr=cross(t1, t2);
if vr(3)<0
    v=v3;
    v3=v1;
    v1=v;
end

% Get the angles of the three vectors
a1=plangle(v1); 
a3=plangle(v3);
if a3<a1
    a3=a3+2*pi;
end

crv=nrlcirc(R, [center,0], a1, a3);    


%% Demo
% pnt1=[0.5,0]; pnt2=[1,1]; pnt3=[0,1]; 
% crv=nrl3pntsarc(pnt1, pnt2, pnt3);
% nrlplot(crv); hold on;
% plot(pnt1(1), pnt1(2), 'ro');
% plot(pnt2(1), pnt2(2), 'ro');
% plot(pnt3(1), pnt3(2), 'ro');





