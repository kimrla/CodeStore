function [d, tt, pntm, h0]=nearestpts(nrl, pts, h0, n, it)

% ANEARESTPTS: Compute the approximated nearest points of a nurve
% 
% Calling Sequence:
% 
%    [d, tt, pntm]=nearestpts(nrl, pts)
%
%    [d, tt, pntm]=nearestpts(nrl, pts, h0)
%
%    [d, tt, pntm]=nearestpts(nrl, pts, [], n, it)
%
%    [d, tt, pntm]=nearestpts(nrl, pts, h0, n, it)
% 
% INPUT:
% 
%    crv	: A NURL curve or surface, see nrlmake.
%
%    pts   : The points.
%
%    h0   : Initial edge length of elements.
%
%    L : Length of the curve or edge of the surface 
%it:迭代次数    n:每次迭代的迭代节点个数
% 
% OUTPUT:
%
%     d - the nearest distance
%
%     tt - parametric points on the curve or surface
%
%     pntm - x, y, z coordinates on the curve or surface
%

if nargin<=3
    n=5;
end

if nargin<=4
    it=5;
end

if (iscell(nrl.knots))
    if (size(nrl.knots,2) == 3)
        %% NURBS structure represents a volume
        
        error('Volume is not surpported!');
        
    elseif (size(nrl.knots,2) == 2)
        %% NURBS structure represents a surface
        
        % Lngth of curves
        if nargin==2 || isempty(h0)
            [M, N, h0]=nrlsrfseednum(nrl);
        elseif nargin>2
            [M, N]=nrlsrfseednum(nrl, h0);
        end
        s=linspace(0, 1, M);
        t=linspace(0, 1, N);

        % Get approximated points
        ds=s(2); dt=t(2);
        [d, tt, pntm]=anearestptsrf(nrl, pts, {s, t});
        k=length(d);
        for j=1:it
            for i=1:k
                s=linspace(max([0, tt(i,1)-ds]), min([1, tt(i,1)+ds]), n);
                t=linspace(max([0, tt(i,2)-dt]), min([1, tt(i,2)+dt]), n);
                [d(i), tt(i,:), pntm(:,i)]=anearestptsrf(nrl, pts(:,i), {s, t});
            end
            ds=ds/2; dt=dt/2;
        end
    end
else
    %% NURBS structure represents a curve

    % Lngth of out circle
    if nargin==2 || isempty(h0)
        [L, h0]=crvmeasure(nrl);
    else
        L = nrlmeasure (nrl);
    end

    % Get the nerest points
    N=round(L/h0);
    s=linspace(0, 1, N); 
    [d, tt, pntm]=anearestptcrv(nrl, pts, s);
    k=length(d); ds=s(2); 
    for j=1:it
        for i=1:k
            s=linspace(max([0, tt(i)-ds]), min([1, tt(i)+ds]), n);
            [d(i), tt(i), pntm(:,i)]=anearestptcrv(nrl, pts(:,i), s);
        end
        ds=ds/2;
    end    
end


%% demo - curve
% % Give a point or a set of points
% pts=[1,1,0; -1,1.1,0; -1.2,-0.6,0; 1,-1,0];
% 
% % Create a polygon
% crv=nrlpolygon(5, 1, [0,0],0); 
% if size(crv,2)>1,
%     crv=nrlglues(crv);
% end
% 
% nrl=crv; pts=pts';
% % Get approximated nerest points
% [d, tt, pntm]=nearestpts(nrl, pts, [], 5, 5);
% 
% % Plot results
% figure; hold on;
% nrlplot(crv, 101);
% plot(pts(1,:), pts(2,:), 'ro');
% plot(pntm(1,:), pntm(2,:), 'r*');


%% demo - curve
% % Give a point or a set of points
% pts=[1,1,0; -1,1.1,0; -1.2,-0.6,0; 1,-1,0];
% 
% % Create a circle
% crv=nrlcirc(1, [0,0], 0, 2*pi);
% 
% nrl=crv; pts=pts';
% % Get approximated nerest points
% [d, tt, pntm]=nearestpts(nrl, pts, [], 5, 5);
% 
% % Plot results
% figure; hold on;
% nrlplot(crv, 101);
% plot(pts(1,:), pts(2,:), 'ro');
% plot(pntm(1,:), pntm(2,:), 'r*');


%% demo - surface
% % Create a surface
% circ=nrlcirc(4, [5,5,5], 0, pi);
% srf=nrlrevolve(circ, [5,5,5], [1,0,0]);
% pts=nrleval(srf, {[0.2, 0.4, 0.8]; [0.2, 0.4, 0.8]});
% pts=pts(:,:);
% 
% figure; 
% nrlplot(srf, [30,30]);
% hold on;
% plot3(pts(1,:), pts(2,:), pts(3,:), 'ro');
% shading interp;
% 
% % Get approximated points
% [d, tt, pntm]=nearestpts(srf, pts, 0.5, 8, 8);
% 
% plot3(pntm(1,:), pntm(2,:), pntm(3,:), 'b*');






