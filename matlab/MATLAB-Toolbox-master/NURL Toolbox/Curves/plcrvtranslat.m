function tcrv=plcrvtranslat(crv, vec)

% Apply a translation to a plane curve
% 
%  Inputs: 
%
%     crv - the curve to be transformed
%     vec - a plane vector like [1.0 0.2]
% 
%  Output: 
%
%      crv - the new curve get by translation
% 
%  Examples:  
% 
%     pnt1=[0, 2]; pnt2=[-1.0, -1.3];
%     crv=nrlline(pnt1, pnt2);
%     vec=[1.5, 1.0];
%     tcrv=plcrvtranslat(crv, vec);
%

tvec = [vec, 0];
st = vectrans(tvec);
tcrv = nrltform(crv, st);


%% Demo
% pnt1=[0, 2]; pnt2=[-1.0, -1.3];
% crv=nrlline(pnt1, pnt2);
% vec=[1.5, 1.0];
% 
% tcrv=plcrvtranslat(crv, vec);
% 
% figure; hold on;
% nrlplot(crv);
% nrlplot(tcrv);


