function crv=nrlquadr(pnts)

% Create a nurl quadrangle formed by 4 points
% 
%  Inputs : pnts - points
%
%  Example   : 
%      pnts={[0.0 0.0 0.3], [1.0 0.0 -0.3], [0.0 1.0 -0.3], [1.0 1.0 0.3]};
%      crv=nrlquadr(pnts);
% 

pntt=sort4nodes(pnts);
crv(1)=nrlline(pntt{1}, pntt{2});
crv(2)=nrlline(pntt{3}, pntt{4});
crv(3)=nrlline(pntt{1}, pntt{3});
crv(4)=nrlline(pntt{2}, pntt{4});


%% Demo
% pnts={[0.0 0.0 0.3], [1.0 0.0 -0.3], [0.0 1.0 -0.3], [1.0 1.0 0.3]};
% crv=nrlquadr(pnts);
% nrlplot(crv(1)); hold on;
% nrlplot(crv(2));
% nrlplot(crv(3));
% nrlplot(crv(4));



