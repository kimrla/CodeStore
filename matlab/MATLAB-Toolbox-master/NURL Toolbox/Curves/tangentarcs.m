function [crvs1, center1, crvs2, center2, crvs]=tangentarcs(lin1, lin2)

% Get the tangent arcs of two straight lines
% 
%  Inputs: 
%
%     crv1, crv2 - the two curves to be used
% 
%  Output: 
%
%     crvs1 - the curves of the small arc and the lines left
%     center1 - center of the small arc
%     crvs2 - the curves of the large arc and the lines left
%     center2 - center of the large arc
%     crvs - the curves of both arcs and the lines left
% 
%  Examples: 
% 
%     lin1=nrlline([0,0], [1,0]);
%     lin2=nrlline([0,0.6], [1,1.1]);
%     [crvs1, center1, crvs2, center2, crvs]=tangentarcs(lin1, lin2);
% 

[crvs1, center1]=tangentarc(lin1, lin2, 1);
[crvs2, center2]=tangentarc(lin1, lin2, 2);   
crvs=tangentarc(crvs1(1), crvs1(2), 2);
n=length(crvs1);
if n==3
    crvs=[crvs, crvs1(3)];
end

%% Demo
% lin1=nrlline([0,0], [1,0]);
% lin2=nrlline([0,0.6], [1,1.1]);
% 
% nrlplot(lin1);
% hold on;
% nrlplot(lin2);
% 
% [crvs1, center1, crvs2, center2, crvs]=tangentarcs(lin1, lin2);
% 
% figure; hold on;
% for i=1:3
%     nrlplot(crvs1(i));
% end
% plot(center1(1), center1(2), 'ro');
% 
% figure; hold on;
% for i=1:3
%     nrlplot(crvs2(i));
% end
% plot(center2(1), center2(2), 'ro');
% 
% figure; hold on;
% for i=1:4
%     nrlplot(crvs(i));
% end
% plot(center1(1), center1(2), 'ro');
% plot(center2(1), center2(2), 'ro');


