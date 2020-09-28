function  crvn = plcrvmirror(crv, line)

% Get the mirror of plane curve(s)
% 
%  Inputs: 
%
%      crv - the curve(s) to be mirrored
% 
%      line - the line to mirror the curve
% 
%  Output: 
%
%      crvn - the new curve(s) get by mirror
% 
%  Examples:  
%     pnt1=[0.8, 0]; pnt2=[1.0, 1.3];
%     crv=nrlline(pnt1, pnt2);
%     line=nrlline([-0.5, 0.5], [1, 2]);
%     tcrv = plcrvmirror(crv, line);
% 

nc=numel(crv);

if nc>1
    crvn(nc) =  plcrvmirror(crv(nc), line); 
    for i=1:nc-1
        crvn(i) =  plcrvmirror(crv(i), line); 
    end
else
    % Check the line
    isl=isnrlline(line);
    if ~isl
        pnts = nrleval(line, [0 1])';
        line=nrlline(pnts(1,1:2), pnts(2,1:2));
    end

    % Get direction of the line
    pnts = nrleval(line, [0 1])';
    nrm=pnts(1,1:2)-pnts(2,1:2);

    % Mirror around the line and the origin or the curve
    vector = [nrm, 0];
    pnt = nrleval(crv, 0)';
    st = vectrans(-pnt);
    crv = nrltform(crv, st);
    rt = vecrot (pi, vector);
    crv = nrltform(crv, rt);
    st = vectrans(pnt);
    crv = nrltform(crv, st);

    % Translate the curve to get its mirror with respect to the line
    ABC= planeline('2points', pnts(1,1:2), pnts(2,1:2));
    d=planeline('distance', pnt(1:2), ABC );
    ABCn= planeline('pointnormal', pnt(1:2), nrm);
    inter=planeline('intersection', ABC, ABCn );
    dr=2*d*vecnorm([inter-pnt(1:2), 0]);
    st = vectrans(dr);
    crv = nrltform(crv, st);
    crvn=crv;
end

%% Demo - line
% pnt1=[0.8, 0]; pnt2=[1.0, 1.3];
% crv=nrlline(pnt1, pnt2);
% line=nrlline([-0.5, 0.5], [1, 2]);
% 
% figure; hold on;
% nrlplot(crv);
% nrlplot(line);
% 
% tcrv = plcrvmirror(crv, line);
% 
% nrlplot(tcrv);


%% Demo - curves
% pnt1=[0.8, 0]; pnt2=[1.0, 1.3];
% crv1=nrlcirc(0.3, [0.6, 0]);
% crv2=nrlcirc(0.5, [0.6, 0]);
% crvs=[crv1, crv2];
% line=nrlline([-0.5, 0.5], [1, 2]);
% 
% figure; hold on;
% nrlplot(crv1, 100);
% nrlplot(crv2, 100);
% nrlplot(line);
% 
% tcrvs = plcrvmirror(crvs, line);
% for i=1:2
%     nrlplot(tcrvs(i), 100);
% end






