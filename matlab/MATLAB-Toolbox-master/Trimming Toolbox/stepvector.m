function sv=stepvector(srf1,srf2,pts1,pts2,pnts,varargin)
% Calculate the step vector in tracing intersection method.

% Input:
%   srf1, srf2: NURBS/Bezier surface structure.
%   pts1,pts2,pnts: Par/phy-coords of the existing precise intersection points,
%       which is a cell array.pts/pnts can be {p} or {p1,p2}, which may be
%       used under the condition that the 2 patches are parallel and the step
%       vector cannot be confirmed by the normal method.
%   'step': Fixed par-domain step length, default is 1e-3.
%   'angle': Fixed step angle, which is used in adaptive step, default is
%       2*pi/1e4.
%   'mode': Fix step or adaptive step, default is 'fix' step length.

% Output:
%   sv: Step vector, including step direction and length. Relative vector.

% Notice: 1.The output sv cannot ensure the correct direction of the step vector,
% maybe the result is the opposite director, in this case just use the
% negative vector of sv, which is -sv. 
% 2. The fixed step is determined in parameter domain instead of physical
% domain, which is helpful to reflect the local geometric information.

STEP=1E-3;
ANGLE=2*pi/1e4;
MODE='fix';
TOL=1E-10;
options=struct('step',STEP,...
    'angle',ANGLE,...
    'tolerance',TOL,...
    'mode',MODE);
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

% Extract the final intersection point in the intersection list.
nump=numel(pnts);
% pts=[pts{:}]; % size(pts)=[2,nump]
pnts=[pnts{:}]; % size(pnts)=[3,nump]
pts1=pts1{end};
pts2=pts2{end};
pnts1=pnts(:,end);
% Calculate the step direction
n1=unitnvector(srf1,pts1);
n2=unitnvector(srf2,pts2);
t=cross(n1,n2);
% The 2 normal vectors are parallel
if norm(t)<options.tolerance
    if nump==2
        t=pnts1-pnts(:,1);
    else % Radial step
        temn=(n1+n2)/2;
        temn=[temn(2),-temn(1),0]; % Relative vector
        temn=temn+pnts1; % Absolute vector
        for j=1:8
            [~,tempnts,d,~]=ite4par(srf1,srf2,[pts1(:),pts2(:)],temn,...
                 'iterative',options.iterative,...
                'mountain',options.mountain,...
                'tolerance',options.tolerance);
            if d<options.tolerance
                t=tempnts-pnts1;
                break;               
            else
                theta=j*2*pi/8;
                rx=vecrot3(pnts1,(n1+n2)/2,theta); 
                temn=rx*[temn;1]; % Absolute vector.
                temn(end)=[];% Homogeneous->Nonhomonegeous
            end    
        end
        if j==8
            sv=[];
            warning('There''s no further step vector or intersection point');
            return;
        end
    end  
end

% Calculate the step length.
if strcmp(options.mode,'fix')
    delta_p=ones(size(pts1))*options.step;
    newpts1=pts1+delta_p;
    ll=newpts1>1;
    newpts1(ll)=pts1(ll)-delta_p(ll);
    newpts2=pts2+delta_p;
    ll=newpts2>1;
    newpts2(ll)=pts2(ll)-delta_p(ll);
    c0=(nrom(nrbeval(srf1,newpts1)-pnts1)+...
        norm(nrbeval(srf2,newpts2)-pnts1))/2;
else % Adaptive step. Use approximate local curvature to determine the step.
    [dsrf11,dsrf12]=nrbderiv(srf1);
    lg1=isderiv(dsrf11);
    lg2=isderiv(dsrf12);
    if lg1==0 && lg2==0
        [~,dF,dF2]=nrbdeval(srf1,dsrf11,dsrf12);
        % 2nd surface invariance II
        L=dot(dF2{1},n1);
        M=dot(dF2{2},n1);
        N=dot(dF2{4},n1);
        % 1st surface invariance I
        E=dot(dF{1},dF{1});
        F=dot(dF{1},dF{2});
        G=dot(dF{2},dF{2});
        K=(L*N-M*M)/(E*F-G*G); % Gauss curvature
        c0=1/K*options.angle;
    else
        
    end
end

sv=c0*vecnorm(t);

end