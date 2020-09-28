function [d, tt, pntm]=anearestpts(nrl, pts, h0, L)

% ANEARESTPTS: Compute the approximated nearest points of a nurve
% 
% Calling Sequence:
% 
%  [d, tt, pntm]=anearestpts(nrl, pts, h0)
%  [d, tt, pntm]=anearestpts(nrl, pts, h0, L)
% 
% INPUT:
% 
%   crv	: A NURL curve or surface, see nrbmak.
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


if (iscell(nrl.knots))
    if (size(nrl.knots,2) == 3)
        %% NURL structure represents a volume
        
        error('Volume is not surpported!');
        
    elseif (size(nrl.knots,2) == 2)
        %% NURL structure represents a surface
        
        % Lngth of curves
        if nargin==3
            crvs=nrlextract(nrl);
            L1 = nrlmeasure (crvs(1));
            L2 = nrlmeasure (crvs(2));
            L3 = nrlmeasure (crvs(3));
            L4 = nrlmeasure (crvs(4));
            Ly=max([L1,L2]);
            Lx=max([L3,L4]);
        elseif nargin==4
            Lx=L; Ly=L;
        end
        
        % Get the nerest points
        M=round(Lx/h0);
        N=round(Ly/h0)+1;
        if M==0 && N>0
            M=N;
        elseif N==0 && M>0
            N=M;
        end
        s=linspace(0, 1, M);
        t=linspace(0, 1, N);
        pnts=nrleval(nrl, {s,t});
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

        % Optimize the nerest points
        n=length(tt);

        for i=1:n
            for j=1:3
                uvi=distpntsrf(nrl, pts(:,i)', tt(i,1), tt(i,2), 1);
                tt(i,:)=uvi(:);                
                if tt(i,1)<0
                    p1=nrleval(nrl, [1+tt(i,1); tt(i,2)]);
                    d1=norm(pts(:,i)-p1);
                    p2=nrleval(nrl, [0; tt(i,2)]);
                    d2=norm(pts(:,i)-p2);
                    if d1<d2
                        tt(i,1)=1+tt(i,1);
                        pntm(:,i)=p1(:);
                        d(i)=d1;
                    else
                        tt(i,1)=0;
                        pntm(:,i)=p2(:);
                        d(i)=d2;
                    end
                elseif tt(i,1)>1
                    p1=nrleval(nrl, [tt(i,1)-1; tt(i,2)]);
                    d1=norm(pts(:,i)-p1);
                    p2=nrleval(nrl, [1; tt(i,2)]);
                    d2=norm(pts(:,i)-p2);
                    if d1<d2
                        tt(i,1)=tt(i,1)-1;
                    else
                        tt(i,1)=1;
                    end
                end
                if tt(i,2)<0
                    p1=nrleval(nrl, [tt(i,1); tt(i,2)+1]);
                    d1=norm(pts(:,i)-p1);
                    p2=nrleval(nrl, [tt(i,1); 0]);
                    d2=norm(pts(:,i)-p2);
                    if d1<d2
                        tt(i,2)=1+tt(i,2);
                    else
                        tt(i,2)=0;
                    end
                elseif tt(i,2)>1
                    p1=nrleval(nrl, [tt(i,1); tt(i,2)-1]);
                    d1=norm(pts(:,i)-p1);
                    p2=nrleval(nrl, [tt(i,1); 1]);
                    d2=norm(pts(:,i)-p2);
                    if d1<d2
                        tt(i,2)=tt(i,2)-1;
                    else
                        tt(i,2)=1;
                    end
                end
                pntm(:,i)=nrleval(nrl, tt(i,:)');
                d(i)=norm(pntm(:,i)-pts(:,i));
            end
        end
    end
else
    %% NURL structure represents a curve

    % Lngth of out circle
    if nargin==3
        L = nrbmeasure (nrl);
    end

    % Get the nerest points
    N=round(L/h0);
    s=linspace(0, 1, N);
    pnts=nrleval(nrl, s);
    dm = DistanceMatrix(pnts', pts');
    [d, id]=min(dm);
    tt=s(id);
    pntm=pnts(:,id);

    % Optimize the nerest points
    n=length(tt);
    for i=1:n
        for j=1:3
            [pi, dpi]=nrldeval(nrl, tt(i));
            dr=pi(:)-pts(:,i);
            tt(i)=tt(i)-dot(dr, dpi)/dot(dpi, dpi);
            if tt(i)<0
                p1=nrleval(nrl, 1+tt(i));
                d1=norm(pts(:,i)-p1);
                p2=nrleval(nrl, 0);
                d2=norm(pts(:,i)-p2);
                if d1<d2
                    tt(i)=1+tt(i);
                    pntm(:,i)=p1(:);
                    d(i)=d1;
                else
                    tt(i)=0;
                    pntm(:,i)=p2(:);
                    d(i)=d2;
                end
            elseif tt(i)>1
                p1=nrleval(nrl, tt(i)-1);
                d1=norm(pts(:,i)-p1);
                p2=nrleval(nrl, 1);
                d2=norm(pts(:,i)-p2);
                if d1<d2
                    tt(i)=tt(i)-1;
                    pntm(:,i)=p1(:);
                    d(i)=d1;
                else
                    tt(i)=1;
                    pntm(:,i)=p2(:);
                    d(i)=d2;
                end
            end
        end
    end
end



%% demo - curve
% % Give a point or a set of points
% pout=[1,1,0; -1,1.1,0; -1,-1,0; 1,-1,0];
% 
% % Mesh seed length
% h0=0.1; 
% 
% % Create a out circle
% crv=nrlpolygon(5, 1, [0,0], 0); 
% if size(crv,2)>1,
%     crv=nrlglues(crv); 
% end
% figure; nrlplot(crv, 100);
% hold on;
% plot(pout(:,1), pout(:,2), 'ro');
% 
% % Lngth of out circle
% L = nrlmeasure (crv);
% 
% % Get the nerest point
% [d, um, pntm]=anearestpts(crv, pout', h0, L);
% 
% plot(pntm(1,:), pntm(2,:), 'r*');

%% demo - surface
% % Give a point or a set of points
% pts=[1,1,5; 9,1,9; 9,9,11; 1,9,6; 5,5,8]';
% 
% % Mesh seed length
% h0=2; 
% 
% % Create a nurbs surface
% srf=nrltestsrf;
% % circ=nrlcirc(5, [5,5,0], 0, pi);
% % srf=nrlrevolve(circ, [5,5,0], [1,0,0]);
% figure; nrlplot(srf, [20, 20]);
% hold on; 
% shading interp;
% plot3(pts(1,:), pts(2,:), pts(3,:), 'ro');
% 
% % Get the nerest points
% [d, uvm, pntm]=anearestpts(srf, pts, h0);
% 
% plot3(pntm(1,:), pntm(2,:), pntm(3,:), 'r*');




