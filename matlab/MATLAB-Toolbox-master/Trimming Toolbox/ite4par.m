function [pts,pnts,d,npn]=ite4par(srf1,srf2,pts0,pnts0,varargin)
% 4-par iterative method.

% Input:
%   srf1,srf2: NURBS/Bezier surface structure.
%   pts0: Parametric coords of the initial point. If size(pts0)=[1,4],
%       pts0=[u0,v0,s0,t0]. If size(pts0)=[1,2], pts0=[p0,q0] which is
%       approximated and should be projected into the 2 surfaces firstly.
%   pnts0: Physical coords of the initial point.size(pnts0)=[3,2] or [2,3].
%       If size(pts0)=[1,2], size(pnts0)=[3,1] or [1,3], which should be 
%       projected into the 2 surfaces firstly.
%   ite: Iterative times, the default equals to 3.
%   method: The default method is Newton-Raphson iterative method And if 
%       input 'mountain', then use Mountain method compulsorily.If the 
%       derivative cannot be calculated (eg. irregular surface), then
%       automatically use Secant method. 
%   
% Output:
%   pts: Precise par-coords of the intersection point. pts=[u,s;v,t].
%   pnts: Precise phy-coords of the intersection point. pnts=[x1,y1,z1;x2,
%       y2,z2]. The final precise point can be obtained by (pnts1+pnts2)/2.
%   d: Distance of the 2 precise points in the 2 surfaces, which is used to
%       determine the error.
%   npn: Parameters to record the positions of the 2 surfaces where the
%       normal vectors are parallel.{[U,V],[S,T]}.

% Default value of the optional variance.
ITE=5;
MIU=[];
TOL=1e-10;

options=struct('iterative',ITE,...
    'mountain',MIU,...
    'tolerance',TOL);
optionnames=fieldnames(options);

nargin=length(varargin);
if round(nargin/2)~=nargin/2
    error('The input is error!');
end

for pair=reshape(varargin,2,[])
    if any(strcmp(pair{1},optionnames))
        options.(pair{1})=pair{2};
    end
end

num1=length(pts0);
npn=[];
if num1==2 % Only one approximate point, which ie needed to be projected. 
    % And replace the glued parameter s* or t*.
    pts1=pts0;
    pts2=pts0;
else % size(pts0)=[1,4]
    pts1=pts0(1:2);
    pts2=pts0(3:4);
end

num2=numel(pnts0);
if num2==3% Only one point, needing to projected into the surfaces respectively.
    pnts0=pnts0(:);
    [pts1,pnts1,~]=projpt2surf(srf1,pts1,pnts0,...
                                'iterative',options.iterative,...
                                'mountain',options.mountain,...
                                'tolerance',options.tolerance,...
                                'order',options.order);
    [pts2,pnts2,~]=projpt2surf(srf2,pts2,pnts0,...
                                'iterative',options.iterative,...
                                'mountain',options.mountain,...
                                'tolerance',options.tolerance,...
                                'order',options.order);
else% num2=6
    m=size(pnts0);
    pnts0=shiftdim(pnts0,find(m~=3));
    pnts1=pnts0(:,1);
    pnts2=pnts0(:,2);
    
end

pts=[pts1(:),pts2(:)];
pnts=[pnts1(:),pnts2(:)];
d=norm(pnts1-pnts2);

for i=1:options.iterative % Endding condition 1.
    n1=unitnvector(srf1,pts1);
    n2=unitnvector(srf2,pts2);
    if abs(dot(n1,n2)-1)<options.tolerance
        warning('The 2 normal vectors are parallel');
        pts=[pts1(:),pts2(:)];
        pnts=[pnts1(:),pnts2(:)];
        d=norm(pnts1-pnts2);
        npn={pts1,pts2};
        break;
    else
        C0=1/2*(1/(sqrt(dot(n1,n2))-1));
        delta_p1=pnts2+pnts1;
        delta_p2=pnts2-pnts1;
        if norm(delta_p2)<options.tolerance
            pts=[pts1(:),pts2(:)];
            pnts=[pnts1(:),pnts2(:)];
            d=norm(pnts1-pnts2);
            break;
        else
            P=delta_p1/2+C0*(dot(delta_p2,n2)*(dot(n1,n2)*n1-n2)-...
                dot(delta_p2,n1)*(dot(n1,n2)*n2-n1));
            [pts1,pnts1,~]=projpt2surf(srf1,(pts1+pts2)/2,P,...
                                 'iterative',options.iterative,...
                                'mountain',options.mountain,...
                                'tolerance',options.tolerance,...
                                'order',options.order);
            [pts2,pnts2,~]=projpt2surf(srf2,(pts1+pts2)/2,P,...
                                'iterative',options.iterative,...
                                'mountain',options.mountain,...
                                'tolerance',options.tolerance,...
                                'order',options.order);
        end
    end

end

end






