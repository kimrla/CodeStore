function tcrv=crvrotate(crv, angle, pnt)

% Apply a rotation to a plane curve around a point
%  
%  Inputs: 
%
%     crv - the curve to be transformed
%     angle - the angle to rotate the curve
%     pnt - a plane pnt like [1.0, 0.2]
% 
%  Output:
%
%     crv - the new curve get by rotation
% 
%  Examples:  
%     pnt1=[0.5, 0]; pnt2=[1.0, 1.3];
%     crv=nrlline(pnt1, pnt2);
%     angle = pi/2; pnt = pnt1; 
%     tcrv=crvrotate(crv, angle, pnt);
% 

if nargin==2
    pnt=[0, 0];
end

tvec = [pnt, 0];
st = vectrans(-tvec);
tcrv = nrltform(crv, st);
rz = vecrotz(angle);
tcrv = nrltform(tcrv, rz);
st = vectrans(tvec);
tcrv = nrltform(tcrv, st);



%% Demo
% pnt1=[0.5, 0]; pnt2=[1.0, 1.3];
% crv=nrlline(pnt1, pnt2);
% angle = pi/2; pnt = pnt1; 
% 
% tcrv=crvrotate(crv, angle, pnt);
% 
% figure; hold on;
% nrlplot(crv);
% nrlplot(tcrv);





