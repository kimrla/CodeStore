function [crvs, center]=nrlfillet(lin1, lin2, R)

% Chamfering two straight lines
% 
% Inputs: 
% 
%     crv1, crv2 - the two curves to be used
% 
%     R - the radius of the circle
% 
%  Output: 
%
%     crvs - the resulted three curves center - center of the arc
% 
%  Examples:  
%     R=0.6;
%     pnt1=[0, 2]; pnt2=[-1.0, -1.3]; pnt3=[2, 0];
%     lin1=nrlline(pnt1, pnt2);
%     lin2=nrlline(pnt2, pnt3);
%     [crvs, center]=nrlfillet(lin1, lin2, R); 
%

% Get the coomon points of two lines
[i, j]=cmpoints(lin1, lin2);
if i(1)==2
    lin1 = nrlreverse(lin1);
end
if j(1)==2
    lin2 = nrlreverse(lin2);
end

% Get the curvature of the lines to judge whether they are straight
[~, c1]=curvature(lin1, 0); 
[~, c2]=curvature(lin2, 0); 
dist1 = nrlmeasure(lin1); 
dist2 = nrlmeasure(lin2); 
if max(abs([c1/dist1, c2/dist2]))>1e-3
    error('One or both lines are not straight at the common ends.');
end

% Extract the two end points of the two lines and get 
%     the directions of the two lines at the two points 
[startendps1, directs1]=nrlcrvextract(lin1);
[~, directs2]=nrlcrvextract(lin2);
pnt=startendps1(1:2, 1); 
dr1=directs1(1:2, 1);
dr2=directs2(1:2, 1);
[dr, q, dr1, dr2]=BisectVector(dr1, dr2);

% Get the center and intersection points
if abs(sin(q))>0 
    L=R/sin(q); S=R/tan(q);
else
    error('The two line are the same.');
end
center=pnt+L*dr';
int1=pnt+S*dr1';
int2=pnt+S*dr2';
int3=center-R*dr';
[~, u1]=crvnearpnt(lin1, [int1; 0]);
[~, u2]=crvnearpnt(lin2, [int2; 0]);
[~, crvs(1)]=nrlsplit(lin1, u1);
[~, crvs(2)]=nrlsplit(lin2, u2);
crvs(3)=nrl3pntsarc(int2', int3', int1');


%% Demo
% R=0.6;
% pnt1=[0, 2]; pnt2=[-1.0, -1.3]; pnt3=[2, 0];
% lin1=nrlline(pnt1, pnt2);
% lin2=nrlline(pnt2, pnt3);
% 
% figure; hold on;
% nrlplot(lin1);
% nrlplot(lin2);
% 
% [crvs, center]=nrlfillet(lin1, lin2, R); 
% 
% figure; hold on;
% for i=1:3
%     nrlplot(crvs(i));
% end





