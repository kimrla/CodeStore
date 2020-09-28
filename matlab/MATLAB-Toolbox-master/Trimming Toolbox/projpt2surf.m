function [pts2,pnts2,d]=projpt2surf(srf,pts1,pnts1,varargin)
% Projection from one given point to the surface, or inverse the parametric
% coords of the given physical coorde, using Newton-Mountain method,
% comparing with nrbsrfreverse.m.

% Three endding conditions: iterative, tolerance, mountain.

% Input: 
%   srf: NURBS/Bezier structure surface.
%   pts1, pnts1: Par/phy-coords of the given point.
%   'iterative':  Iimes of iterative, default is 3.
%   'order': 1 for 1-degree derivative (default) and 2 for 2-degree derivative
%       (more precise and fast). Once the derivatie (both 1 and 2-degree),
%       the derivative information is replaced by the difference quotient.
%   'tolerance'=TOL: Tolerance between x(k+1) and x(k), default is eps.
%   'mountain'=MIU: Mountain factor, default is 1/2e10

% Output:
%   pts2: Par-coords of the prjoection point in the surface.
%   pnts2: Phy-coords of the prjoection point in the surface.
%   d: Distance between the 2 points.

% Default value of the optional variance.
ITE=3;
MIU=1/2e10;
TOL=eps;
ORD=1;

options=struct('iterative',ITE,...
    'mountain',MIU,...
    'tolerance',TOL,...
    'order',ORD);
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

if options.order==1  
    dsrf=nrbderiv(srf);
    lg=isderiv(dsrf);    
    if lg==1 % The derivative information needs to be replaced by 
        % difference quotient       
        dr=@(knt) (nrbeval(srf,knt)-pnts1);
        dF=discretederiv(dr,pts1);
        dF=dF'; % size(dF)=[2,3]        
        f=@(knt) (dF*dr(knt));              
        [pts2,pnts2,d]=NewtonDiscrete(f,pts1,'iterative',ITE,...
                                            'mountain',MIU,...
                                            'tolerance',TOL);             
    else    
        % dr=nrbeval(srf,knt)-pnts1), size(dr)=[3,numpoints]
        % Handle of original function, size(f(x))=[2,1], size(df(x))=[2,2]
        f=@(knt) (jacobi(srf,dsrf,knt)*(nrbeval(srf,knt)-pnts1));
        df=@(knt) (jacobi(srf,dsrf,knt)*jacobi(srf,dsrf,knt)');
        [pts2,pnts2,d]=Newton(f,df,pts1,'iterative',ITE,...
                                            'mountain',MIU,...
                                            'tolerance',TOL);     
    end       
elseif options.order==2
    dsrf=nrbderiv(srf);
    lg1=isderiv(dsrf); % Check the 1st derivative structure
    dsrf2=nrbderiv(dsrf);
    lg2=isderiv(dsrf2);% Check the 2nd derivative structure 
    if lg1==1 % Both 1st and 2nd derivative are NOT calculated
        dr=@(knt) (nrbeval(srf,knt)-pnts1);
        [dF,dF2]=discretederiv(dr,pts1); % size(dF)=[3,2], size(dF2)=[3,4]
        dF=dF';% size(dF)=[2,3]           
        f=@(knt) (dF*dr(knt));           
        dF=dF.*dF'+[dot(dr(pts1),dF2(:,1)),dot(dr(pts1),dF2(:,2));...
                    dot(dr(pts1),dF2(:,3)),dot(dr(pts1),dF2(:,4))];
        for i=1:options.iterative % Endding condition 1.
            F=f(x0);
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
    elseif lg2==1 % 1st derivative is OK, 2nd derivative is wrong.    
        % dr=nrbeval(srf,knt)-pnts1, size(dr)=[3,numpoints]
        % Handle of original function, size(f(x))=[2,1], size(df(x))=[2,2]
        dr=@(knt) (nrbeval(srf,knt)-pnts1);
        f=@(knt) (jacobi(srf,dsrf,knt)*dr(knt));
        df=@(knt) (jacobi(srf,dsrf,knt)*jacobi(srf,dsrf,knt)');
        pts2=NewtonDiscrete(df,pts1,'iterative',ITE,...
                                            'mountain',MIU,...
                                            'tolerance',TOL);  
        pnts2=f(pts2);
        d=norm(pnts2-pnts1);
                               
    else % Use 1st and 2nd derivative directly
        dr=@(knt) (nrbeval(srf,knt)-pnts1);
        [dF,dF2]=jacobi(srf,dsrf,knt,dsrf2);
        f=@(knt) (dF*dr(knt));
        df=dF*dF';
        df=df+[dot(dr(knt),dF2(1,:)),dot(dr(knt),dF2(2,:));...
                dot(dr(knt),dF2(3,:)),dot(dr(knt),dF2(4,:))];
        [pts2,pnts2,d]=Newton(f,df,pts1,'iterative',ITE,...
                                    'mountain',MIU,...
                                    'tolerance',TOL);     
    end  
else
    error('The derivative is wrong!');
end

end

function [dF,dF2]=jacobi(srf,dsrf,knt,varargin)
if nargout==1
    [~,dF]=nrbdeval(srf,dsrf,knt);
elseif nargout==2
    dsrf2=varargin{2};
    [~,dF,dF2]=nrbdeval(srf,dsrf,knt,dsrf2);
end
if iscell(dF)
    dF=[dF{:}];% size dF=[3,2]
end
if iscell(dF2)
    dF2=[dF2{:}];% size dF=[3,2]
end
dF=dF'; % size(dF)=[2,3] or [1,3]
dF2=dF2';% size(dF2)=[4,3] or [1,3]

end




