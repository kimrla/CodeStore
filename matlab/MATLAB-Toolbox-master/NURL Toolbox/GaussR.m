function [x,C]=GaussR(n,a,b,R)

% Solve the nodes and weights of Gauss integration method
%     
% Calling Sequence:
% 
%   [x, C]=GaussR(knots, order)
%
%  Input:
%    a, b - the span of the interval
%    N - the number of nodes required
%    R - the number of iterations (default is 5)
%
%  Output:
%    x - Gauss nodes
%    C - Gauss weights
%
% Example: 
% 
%   n=12; 
%   [x,C]=GaussR(n,-1,1);
%   f=sqrt(x+1.5); ef=2.399529;
%   intf=sum(C.*f);
%

    if nargin==3
        R=5;    
    end
    if nargin==1
        R=5; 
        a=-1; b=1;
    end

    if n==1
        x=0; C=2;
    elseif n==2
        x=[-0.577350269189626;0.577350269189626];
        C=[1; 1];
    elseif n==3;
        x=[-0.774596669241483;0;0.774596669241483];
        C=[0.555555555555556;0.888888888888889;0.555555555555556];
    elseif n==4
        x=[-0.861136311594053;-0.339981043584856;0.339981043584856;0.861136311594053];
        C=[0.347854845137454;0.652145154862546;0.652145154862546;0.347854845137454];
    elseif n==5
        x=[-0.906179845938664;-0.538469310105683;0;0.538469310105683;0.906179845938664];
        C=[0.236926885056189;0.478628670499367;0.568888888888889;0.478628670499367;0.236926885056189];
    elseif n==6
        x=[-0.932469514203152;-0.661209386466265;-0.238619186083197;0.238619186083197;0.661209386466265;0.932469514203152];
        C=[0.171324492379170;0.360761573048138;0.467913934572691;0.467913934572691;0.360761573048138;0.171324492379170];
    elseif n==7
        x=[-0.949107912342759;-0.741531185599395;-0.405845151377397;0;0.405845151377397;0.741531185599395;0.949107912342759];
        C=[0.129484966168870;0.279705391489277;0.381830050505119;0.417959183673469;0.381830050505119;0.279705391489277;0.129484966168870];
    elseif n==8
        x=[-0.960289856497536;-0.796666477413627;-0.525532409916329;-0.183434642495650;0.183434642495650;0.525532409916329;0.796666477413627;0.960289856497536];
        C=[0.101228536290375;0.222381034453374;0.313706645877887;0.362683783378362;0.362683783378362;0.313706645877887;0.222381034453374;0.101228536290375];
    elseif n==9
        x=[-0.968160239507626;-0.836031107326636;-0.613371432700591;-0.324253423403809;0;0.324253423403809;0.613371432700591;0.836031107326636;0.968160239507626];
        C=[0.0812743883615750;0.180648160694857;0.260610696402935;0.312347077040003;0.330239355001260;0.312347077040003;0.260610696402935;0.180648160694857;0.0812743883615750];
    elseif n==10
        x=[-0.973906528517172;-0.865063366688984;-0.679409568299024;-0.433395394129247;-0.148874338981631;0.148874338981631;0.433395394129247;0.679409568299024;0.865063366688984;0.973906528517172];
        C=[0.0666713443086869;0.149451349150581;0.219086362515982;0.269266719309996;0.295524224714753;0.295524224714753;0.269266719309996;0.219086362515982;0.149451349150581;0.0666713443086869];
    else        
        rLp=LegendreRt(n,R);
        x=rLp;
        vLp=LegendreRecursion(n-1,rLp);
        [~,vdLp]=LegendreRecDer(n,rLp);
        C=(2/n)./(vLp.*vdLp);
    end

    x=(b-a)*x/2+(b+a)/2;
    C=(b-a)*C/2;
end

% Get the roots of Legendre polynomials
%
% rx :    the roots
% rry :  the value of Legendre polynomials at the roots, their should be
%             zeros

function [rx,rry,xr,ry]=LegendreRt(m,R)

    if nargin==1
        R=5;
    end

    n=7*m;
    x=LobattoChebyshev(-1,1,n);
    y=LegendreRecursion(m,x);

    % Get the span of roots
    s=1; 
    ry=zeros(m,2);
    xr=ry; 
    ni=ry(:,1); 
    for i=1:fix(n/2)
        if y(i)*y(i+1)<0
            ry(s,1)=y(i);
            ry(s,2)=y(i+1);
            xr(s,1)=x(i);
            xr(s,2)=x(i+1);
            ni(s)=i;
            s=s+1;
        end
    end

    % Solve the roots of Legendre
    rx=zeros(m,1); rry=rx;
    for i=1:fix(m/2)
        if abs(ry(i,1))<abs(ry(i,2))
            rx(i)=xr(i,1);
        else
            rx(i)=xr(i,2);
        end
        for j=1:R
            [yj,dyj,~]=LegendreRecDer(m,rx(i));
            rx(i)=rx(i)-yj/dyj;
        end
        [rry(i),~,~]=LegendreRecDer(m,rx(i));
    end

    % Expand
    if mod(m,2)==0
        p=fix(m/2);
        for i=1:p
            rx(p+i)=-rx(p-i+1);
            rry(p+i)=rry(p-i+1);
            xr(p+i,1)=-xr(p-i+1,1);
            ry(p+i,1)=ry(p-i+1,1);
            xr(p+i,2)=-xr(p-i+1,2);
            ry(p+i,2)=ry(p-i+1,2);
        end
    elseif mod(m,2)==1
        p=ceil(m/2);
        for i=1:fix(m/2)
            rx(p+i)=-rx(p-i);
            rry(p+i)=rry(p-i);
            xr(p+i,1)=-xr(p-i,1);
            ry(p+i,1)=-ry(p-i,1);
            xr(p+i,2)=-xr(p-i,2);
            ry(p+i,2)=-ry(p-i,2);
        end
    end
end

% Legendre by the recursion formula

function y=LegendreRecursion(n,x)

    if n==0
        y=1;
    elseif n==1
        y=x;
    elseif n>1    
        y1=1; y2=x;
        for k=1:n-1
            y=((2*k+1)/(k+1))*x.*y2-(k/(k+1))*y1;
            y1=y2;
            y2=y;
        end
    end
end

% Get derivatives of Legendre polynomial by the recursion formula

function [y,dy,ddy]=LegendreRecDer(n,x)

    if n==0
        y=1; dy=0; ddy=0;
    elseif n==1
        y=x; dy=1; ddy=0;
    elseif n>1    
        y1=1; y2=x;
        dy1=1;
        ddy1=0;
        for k=1:n-1
            y=((2*k+1)/(k+1))*x.*y2-(k/(k+1))*y1;
            y1=y2;
            y2=y;

            dy=x.*dy1+(k+1)*y1;
            dy0=dy1;
            dy1=dy;

            ddy=x.*ddy1+(k+2)*dy0;
            ddy1=ddy;
        end
    end
end

% Lobatto Chebyshev nodes
function x=LobattoChebyshev(a,b,N)

    x=zeros(N,1); % 求节点坐标
    for j=1:N
        nd=(j-1)/(N-1); h=0.5*(1-cos(nd*pi));
        x(j)=a+h*(b-a);
    end
end







