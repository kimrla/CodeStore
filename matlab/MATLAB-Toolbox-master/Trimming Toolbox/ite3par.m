function [pts,pnts,d]=ite3par(srf1,srf2,st2,pts0,pnts0,varargin)
% 3-par iterative method. The 2 surfaces are S1(u,v) and S2(s,t), where
% parameter s or t is glued (depending on the input) . If the derivative 
% cannot be calculated, then use the difference quotient automatically.

% Input:
%   srf1,srf2: NURBS/Bezier surface structure.
%   st2: [Signal, s or t] of parameter s or t, to determine which parameter 
%       is glued. The first elemtn 1 for gluing the 1st par(s), 2 for 
%       gluing the 2nd par(t). The 2nd element is the glued parameter.
%   pts0: Parametric coords of the initial point. If size(pts0)=[1,3],
%       pts0=[u0,v0,s0 or t0]. If size(pts0)=[1,2], pts0=[p0,q0] which is
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
%       y2,z2].
%   d: Distance of the 2 precise points in the 2 surfaces, which is used to
%       determine the error.

% Default value of the optional variance.
ITE=3;
MIU=[];
TOL=eps;

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
if num1==2 % Only one approximate point, which ie needed to be projected. 
    % And replace the glued parameter s* or t*.
    pts1=pts0;
    pts2=pts0;
else % size(pts0)=[1,3]
    pts1=pts0(1:2);
    pts2=pts0(3);
end
% Set par-coords of srf2. 
if st2(1)==1 % s is glued, s*
    pts2=[st2(2),pts2];    
    st2(1)=2;   
elseif st2(2)==2 % t is glued, t*
    pts2=[pts2,st2(2)];    
    st2(1)=1; 
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


% pts=[pts1(:);pts2(:)];
% pts=pts'; % size(pts)=[1,4]
f=@(pts0) (nrbeval(srf1,pts0(1:2))-nrbeval(srf2,pts2));

if detercond('tolerance',pnts1,pnts2,options.tolerance)
    pts=[pts1(:),pts2(:)];
    pnts=[pnts1(:),pnts2(:)];
    d=norm(pnts1-pnts2);  
    return;
else
    dsrf1=nrbderiv(srf1);
    lg1=isderiv(dsrf1);   
    dsrf2=nrbderiv(srf2);
    lg2=isderiv(dsrf2);  
    if lg1==1 || lg2==1% The derivative information needs to be replaced by 
        % difference quotient   
        x0=pts0(:);
        x0=x0';
        for i=1:options.iterative % Endding condition 1.
            dF1=discretederiv(@(knt) (nrbeval(srf1,knt)),x0(1:2));
            if st2(1)==1 % s is glued, s*
                dF2=discretederiv(@(knt) (nrbeval(srf2,knt)),[x0(3),st2(2)]);
            elseif st2(2)==2 % t is glued, t*
                dF2=discretederiv(@(knt) (nrbeval(srf2,knt)),[st2(2),x0(3)]);
            end    
            dF=[dF1,dF2(:,st2(1))];                       
            % Ending condition 2.    
            if ~isempty(options.mountain)
                x=itemount(f,x0,F,dF,options.mountain);   
            else
                x=x0-dF\F; 
            end     
            % Ending conditionn 3.
            if detercond('tolerance',x0,x,options.tolerance)
                break;   
            else
                x0=x;
            end        
        end   
        if st2(1)==1 % s is glued, s*
            pts=[x(1),st2(2);
                x(2),x(3)];
        elseif st2(2)==2 % t is glued, t*
            pts=[x(1),x(3);
                x(2),st2(2)];
        end
        pnts(:,1)=nrbeval(srf1,pts(:,1));
        pnts(:,2)=nrbeval(srf2,pts(:,2));      
    else % lg1==0,lg2==0    
        % dr=nrbeval(srf,knt)-pnts1), size(dr)=[3,numpoints]
        % Handle of original function, size(f(x))=[2,1], size(df(x))=[2,2]      
        tem4=jacobi(srf2,dsrf2,pts2);
        df=@(pts0) ([jacobi(srf1,dsrf1,pts0(1:2)),tem4(:,st2(1))]); % [Ru,Rv,Rs 
        % or Rt], size(df)=[3,3]
        % 代入牛顿法中，每次传递的是函数（要达到的效果）还是一个不变的值？
        [pts2,~,d]=Newton(f,df,pts0,'iterative',ITE,...
                                            'mountain',MIU,...
                                            'tolerance',TOL);  
        if st2(1)==1 % s is glued, s*
            pts=[pts2(1),st2(2);
                pts2(2),pts2(3)];
        elseif st2(2)==2 % t is glued, t*
            pts=[pts2(1),pts2(3);
                pts2(2),st2(2)];
        end
        pnts(:,1)=nrbeval(srf1,pts(:,1));
        pnts(:,2)=nrbeval(srf2,pts(:,2));                                                             
    end       

end   
pts=pts';
pnts=pnts';% size(pnts)=[2,3]
end

function dF=jacobi(srf,dsrf,knt)

[~,dF]=nrbdeval(srf,dsrf,knt);
if iscell(dF)
    dF=[dF{:}];% size dF=[3,2]
end

end
