function [d, tt, pntm]=nearestpts(nrb, pts, h0, L)

% ANEARESTPTS: Compute the approximated nearest points of a nurve
% 
% Calling Sequence:
% 
%  [d, tt, pntm]=anearestpts(nrb, pts, h0)
%  [d, tt, pntm]=anearestpts(nrb, pts, h0, L)
% 
% INPUT:
% 
%   crv	: A NURBS curve or surface, see nrbmak.
%   pts   : The points.
%   h0   : Initial edge length of elements.
%   L : Length of the curve or edge of the surface 
% 
% OUTPUT:
%
%    d - the nearest distance
%    tt - parametric points on the curve or surface
%    pntm - x, y, z coordinates on the curve or surface
%

if (iscell(nrb.knots))
    if (size(nrb.knots,2) == 3)
        %% NURBS structure represents a volume
        
        error('Volume is not surpported!');
        
    elseif (size(nrb.knots,2) == 2)
        %% NURBS structure represents a surface
        
        % Lngth of curves
        if nargin==3
            crvs=nrbextract(nrb);
            L1 = nrbmeasure (crvs(1));
            L2 = nrbmeasure (crvs(2));
            L3 = nrbmeasure (crvs(3));
            L4 = nrbmeasure (crvs(4));
            Ly=max([L1,L2]);
            Lx=max([L3,L4]);
        elseif nargin==4
            Lx=L; Ly=L;
        end
        
        % Get nodes
        M=round(Lx/h0);
        N=round(Ly/h0)+1;
        if M==0 && N>0
            M=N;
        end
        if N==0 && M>0
            N=M;
        end
        N=N+1;
        s=linspace(0, 1, M);
        t=linspace(0, 1, N);

        % Get approximated points
        ds=s(2); dt=t(2);
        [d, tt, pntm]=nearestptsrf(nrb, pts, {s, t});
        k=length(d); n=6;
        for j=1:3
            for i=1:k
                s=linspace(max([0, tt(i,1)-ds]), min([1, tt(i,1)+ds]), n);
                t=linspace(max([0, tt(i,2)-dt]), min([1, tt(i,2)+dt]), n);
                [d(i), tt(i,:), pntm(:,i)]=nearestptsrf(nrb, pts(:,i), {s, t});
            end
            ds=ds/(n-1); dt=dt/(n-1);
        end

    end
else
    %% NURBS structure represents a curve

    % Lngth of out circle
    if nargin==3
        L = nrbmeasure (nrb);
    end

    % Get the nerest points
    N=round(L/h0);
    s=linspace(0, 1, N); 
    [d, tt, pntm]=nearestptcrv(nrb, pts, s);
    k=length(d); ds=s(2); n=6;
    for j=1:3
        for i=1:k
            s=linspace(max([0, tt(i)-ds]), min([1, tt(i)+ds]), n);
            [d(i), tt(i), pntm(:,i)]=nearestptcrv(nrb, pts(:,i), s);
        end
        ds=ds/(n-1);
    end    
end


function [d, s, pntm]=nearestptcrv(crv, pts, s)

pnts=nrbeval(crv, s);
dm = DistanceMatrix(pnts', pts');
[d, id]=min(dm);
s=s(id);
pntm=pnts(:,id);

function [d, tt, pntm]=nearestptsrf(srf, pts, st)

s=st{1}; t=st{2}; M=length(s);
pnts=nrbeval(srf, st);
dm = DistanceMatrix(pnts(:,:)', pts');
[d, id]=min(dm);
k=length(d);
tt=zeros(k,2);
for i=1:k
    tt(i,1)=rem(id(i),M);
    tt(i,2)=fix(id(i)/M)+1;
end
pp=tt==0;
tt(pp)=M;
tt(pp(:,1),2)=tt(pp(:,1),2)-1;
for i=1:k
    tt(i,1)=s(tt(i,1));
    tt(i,2)=t(tt(i,2));
end
pntm=pnts(:,id);

%% demo - curve
% % Give a point or a set of points
% pout=[1,1,0; -1,1.1,0; -1,-1,0; 1,-1,0];
% 
% % Mesh seed length
% h0=0.1; 
% 
% % Create a out circle
% crv=CreatCurves('polygon', 5, 1, [0,0],0); 
% if size(crv,2)>1,
%     crv=CurveManage('glues',crv);
% end
% figure; nrbplot(crv, 100);
% hold on;
% plot(pout(:,1), pout(:,2), 'ro');
% 
% % Lngth of out circle
% L = nrbmeasure (crv);
% 
% % Get the nerest point
% [d, tt, pntm]=nearestpts(crv, pout', h0);


%% demo - surface
% % Give a point or a set of points
% pts=[1,1,5; 9,1,9; 9,9,9; 1,9,6; 5,5,10]';
% 
% % Mesh seed length
% h0=0.5; 
% 
% % Create a surface
% circ=nrbcirc(4, [5,5,5], 0, pi);
% srf=nrbrevolve(circ, [5,5,5], [1,0,0]);
% 
% figure; 
% nrbplot(srf, [30,30]);
% hold on;
% plot3(pts(1,:), pts(2,:), pts(3,:), 'ro');
% shading interp;
% 
% % Get approximated points
% [d, tt, pntm]=nearestpts(srf, pts, h0);
% 
% plot3(pntm(1,:), pntm(2,:), pntm(3,:), 'r*');






