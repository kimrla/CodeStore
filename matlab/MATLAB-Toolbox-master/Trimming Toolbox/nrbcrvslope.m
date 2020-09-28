function [ pts,pnts,tpnts ] = nrbcrvslope( nrb )
% Get the special points of NURBS curve whose slopes' coordinates are
% x(t)=0 or y(t)=0.

% Input:
%     nrb: NURBS curve structure.
%  Output:
%     pts: Parameter coordinates of special points.
%     pnts: Physical coordinates of special points.
%     tpnts: Tangent vector of special points.

[dnrb,ddnrb]=nrbderiv(nrb);
pts=[];
pnts=[];
tpnts=[];
tol=1e-6; % Tolerance.
iter=7; % Iterative times.

t=linspace(0,1,100);
[~,tnrb]=nrbdeval(nrb,dnrb,t);
tnrb=vecnorm(tnrb);

testnrb=tnrb(:,2:end).*tnrb(:,1:end-1);
for i=1:2
    testid{i}=find(testnrb(i,:)<0);
end

% testx=x(2:end).*x(1:end-1);
% testid=find(testx<0);
% 
% [~,idxy{1}]=sort(abs(x));
% [~,idxy{2}]=sort(abs(y));
% for i=1:2
%     txy{i}=t(idxy{i});
% end

for i=1:2 % x and y cords.
    num=length(testid{i});
    for j=1:num
        it=1;
        id=testid{i}(j);
        tk=(t(id)+t(id+1))/2; % Initial parameter value of Newton method.
%         det1=false;
        while it<=iter % At most iterations.
            [~,f,df]=nrbdeval(nrb,dnrb,ddnrb,tk);
            tk1=tk-f(i)/df(i);
            tk1=min(tk1,1);tk1=max(tk1,0);
            if (abs(tk1-tk)<tol)
                if (abs(f(i)/df(i))<tol) % Ratio of slope to curvature.
                    pts=[pts,tk1]; % Only the special slope points, except the stagnation ponts.
                end
%                 det1=true;          
                break;
            else
                tk=tk1;
            end
        end
%         if (~det1 && det2>1)
%         % Determine each knot sorted by the tangent vector. If pts does no
%         % longer increase, which means the rest larger knots are
%         % impossibily the special points, the break the cycle.
%             break;
%         end
    end
end

if ~isempty(pts) % Remove par cords which are too similar with each other.
    pts=unique(pts);
    tempts=pts(2:end)-pts(1:end-1);
    id=find(tempts<tol);
    pts(id)=[];
    pnts=nrbeval(nrb,pts);
    dnrb=nrbderiv(nrb);
    [~,tpnts]=nrbdeval(nrb,dnrb,pts);
    tpnts=vecnorm(tpnts);
end
    
end

