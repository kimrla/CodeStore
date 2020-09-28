function [pnts, cvt, tors, dp, ddp]=nrlcvttors(crv, t)

% Get the curvature and torsion of a nurl curve
% 
%  Input: 
% 
%    crv - a nurl curve
%    t - the points of evaluation
%
%  Output:
% 
%    pnts - points of the curve
%    cvt - curvature of the curve
%    tors - torsion of the curve
%    dp - first order derivatives of the curve
%    ddp - second order derivatives of the curve
% 

t=checktt(t);
dps = nrldseval (crv, t, 3); 
r12=veccross(dps{2}, dps{3}); 
dr=zeros(size(t)); 
dr12=dr; 
dr123=vecdot(r12, dps{4});
n=length(t); 
for i=1:n 
    dr(i)=norm(dps{2}(:,i)); 
    dr12(i)=norm(r12(:,i)); 
end
cvt=dr12./dr.^3;
tors=dr123./dr12.^2;
pnts=dps{1};
pp=isnan(tors); tors(pp)=0;
dp=dps{2};
ddp=dps{3};

%% Demo 1
% % Constants of the cylindrical spiral curve
% a=1; b=0.1; t1=0; t2=5*pi;
% 
% % Get the curve
% crv=spiralcylind(a, b, t1, t2, 10);
% nrlctrlplot(crv);
% 
% % Get the curvature and torsion of the curve
% t=linspace(0, 1, crv.number+5); 
% [pnts, cvt, tors]=nrlcvttors(crv, t);


%% Demo 2
% % Constants of the cylindrical spiral curve
% a=1; b=0.1; t1=0; t2=5*pi;
% 
% % Get the curve
% crv=nrlcirc;
% nrlctrlplot(crv);
% 
% % Get the curvature and torsion of the curve
% t=linspace(0, 1, crv.number+5); 
% [pnts, cvt, tors]=nrlcvttors(crv, t);



