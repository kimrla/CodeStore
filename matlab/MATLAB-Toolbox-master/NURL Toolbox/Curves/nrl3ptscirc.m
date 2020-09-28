function crv=nrl3ptscirc(pnt1, pnt2, pnt3)

% Create a nurl circle by three points
%
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

center= planeline('center3ps',pnt1,pnt2,pnt3);
R=norm(pnt1-center);
crv=nrlcirc(R, [center,0]);

%% ! Demo
% pnt1=[0.5,0]; pnt2=[1,1]; pnt3=[0,1]; 
% crv=nrl3ptscirc(pnt1, pnt2, pnt3);
% nrlplot(crv);





