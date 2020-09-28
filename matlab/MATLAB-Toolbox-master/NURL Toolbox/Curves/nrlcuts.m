function [crvn, crvsn]=nrlcuts(crv, crvs)

% Cut curves by one curve 
% 
%  Inputs: 
%
%     crv - the curve to cut
% 
%     crvs - the curves to be cutted
% 
%  Output: 
%
%     crvn - the splitted curvs of the curve to cut
% 
%     crvsn - the curves cutted from crvs
% 
%  Examples:  
%     crv=nrlline([3.0, 6], [6, 1.0]);
%     p=2; pnt=[6;5.4;0]; 
%     crv1=nrlcirc(3, pnt', pi, 2*pi);
%     x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
%     y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
%     crv2=nrlspline(p,x,y); 
%     crvs=[crv1, crv2];
%     [crvn, crvsn]=nrlcuts(crv, crvs);
% 

n=numel(crvs);
crvsn=[]; u=[];
for i=1:n
    [~, crvs2, uu]=intsctpnts(crv, crvs(i));
    crvsn=[crvsn, crvs2];
    if ~isempty(uu)  
        u=[u, uu(:,1)'];
    end
end
crvn=nrlsplits(crv, u);


%% Demo
% crv=nrlline([3.0, 6], [6, 1.0]);
% p=2; pnt=[6;5.4;0]; 
% crv1=nrlcirc(3, pnt', pi, 2*pi);
% x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
% y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
% crv2=nrlspline(p,x,y); 
% crvs=[crv1, crv2];
% 
% figure; hold on;
% for i=1:numel(crvs)
%     nrlplot(crvs(i), 100);
% end
% nrlplot(crv, 100);
% 
% [crvn, crvsn]=nrlcuts(crv, crvs);
% 
% figure; hold on;
% for i=1:numel(crvsn)
%     nrlplot(crvsn(i), 100);
% end
% for i=1:numel(crvn)
%     nrlplot(crvn(i), 100);
% end




