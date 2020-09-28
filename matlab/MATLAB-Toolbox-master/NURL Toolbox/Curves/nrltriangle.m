function crv=nrltriangle(pnts)

% Create a nurl triangle formed by three points
%
%  Inputs : 
%     pnts - points
%
%  Example   : 
%     pnts={[0.0 0.0 0.3], [1.0 0.0 -0.3], [0.0 1.0 -0.3]};
%     crv=nrltriangle(pnts);
% 

crv(1)=nrlline(pnts{1}, pnts{2});
crv(2)=nrlline(pnts{2}, pnts{3});
crv(3)=nrlline(pnts{3}, pnts{1});


%% Demo
% pnts={[0.0 0.0 0.3], [1.0 0.0 -0.3], [0.0 1.0 -0.3]};
% crv=nrltriangle(pnts);
% nrlplot(crv(1));
% hold on;
% nrlplot(crv(2));
% nrlplot(crv(3));



