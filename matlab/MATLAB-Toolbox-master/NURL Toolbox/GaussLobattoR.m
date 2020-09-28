function [x,C]=GaussLobattoR(n,a,b,R)
%
% Solve roots and weights of Gauss-Lobatto quadrature
% 
% Calling Sequence:
% 
%   [x,C]=GaussLobattoR(n,a,b)
% 
%   [x,C]=GaussLobattoR(n,a,b,R)
%
% Input : 
% 
%    n - the number of integration nodes
%
%    [a, b] - the interval of integration
%
%    R - the number of interations of solve the nodes
%         Default value is 5
%  
% Output :
% 
%    x - the integration node
%
%    C - the integration weights
% 
% Example: 
% 
%   n=12; 
%   [x,C]=GaussLobattoR(n,-1,1);
%   f=sqrt(x+1.5); ef=2.399529;
%   intf=sum(C.*f);
%

if nargin==3
    R=5;
end
switch n
    case 2
        x=[-1;1];
        C=[1;1];
    case 3
        x=[-1;0;1];
        C=[0.333333333333333;1.33333333333333;0.333333333333333];
    case 4
        x=[-1;-0.447213595499958;0.447213595499958;1];
        C=[0.166666666666667;0.833333333333334;0.833333333333334;0.166666666666667];
    case 5
        x=[-1;-0.654653670707977;0;0.654653670707977;1];
        C=[0.100000000000000;0.544444444444444;0.711111111111111;0.544444444444444;0.100000000000000];
    case 6
        x=[-1;-0.765055323929465;-0.285231516480645;0.285231516480645;0.765055323929465;1];
        C=[0.0666666666666667;0.378474956297847;0.554858377035486;0.554858377035486;0.378474956297847;0.0666666666666667];
    case 7
        x=[-1;-0.830223896278567;-0.468848793470714;-4.07831529249908e-56;0.468848793470714;0.830223896278567;1];
        C=[0.0476190476190476;0.276826047361566;0.431745381209863;0.487619047619048;0.431745381209863;0.276826047361566;0.0476190476190476];
    case 8
        x=[-1;-0.871740148509607;-0.591700181433142;-0.209299217902479;0.209299217902479;0.591700181433142;0.871740148509607;1];
        C=[0.0357142857142857;0.210704227143506;0.341122692483504;0.412458794658704;0.412458794658704;0.341122692483504;0.210704227143506;0.0357142857142857];
    case 9
        x=[-1;-0.899757995411460;-0.677186279510738;-0.363117463826178;0;0.363117463826178;0.677186279510738;0.899757995411460;1];
        C=[0.0277777777777778;0.165495361560805;0.274538712500162;0.346428510973046;0.371519274376417;0.346428510973046;0.274538712500162;0.165495361560805;0.0277777777777778];
    case 10
        x=[-1;-0.919533908166459;-0.738773865105505;-0.477924949810445;-0.165278957666387;0.165278957666387;0.477924949810445;0.738773865105505;0.919533908166459;1];
        C=[0.0222222222222222;0.133305990851070;0.224889342063127;0.292042683679684;0.327539761183898;0.327539761183898;0.292042683679684;0.224889342063127;0.133305990851070;0.0222222222222222];
    otherwise
        rdLp=dLegendreRt(n-1,R);
        vLp=LegendreRecursion(n-1,rdLp);
        x=[-1;rdLp;1];
        C=(2/(n*(n-1)))*[1;1./(vLp.^2);1];
end

x=(b-a)*x/2+(b+a)/2;
C=(b-a)*C/2;

end

% Get the roots of first order derivative of Legendre polynomials
%
% rxd :    the roots of first order derivative of Legendre polynomials
% m :  the order of Legendre polynomials

function rxd=dLegendreRt(m,R)

    if nargin==1
        R=5;
    end

    n=7*m;
    x=LobattoChebyshev(-1,1,n);
    [~,dy,~]=LegendreRecDer(m,x);

    % Get the span of roots
    t=1;
    rdy=zeros(m-1,2); xdr=rdy;
    nid=rdy(:,1);
    for i=1:fix(n/2)
        if dy(i)*dy(i+1)<0
            rdy(t,1)=dy(i);
            rdy(t,2)=dy(i+1);
            xdr(t,1)=x(i);
            xdr(t,2)=x(i+1);
            nid(t)=i;
            t=t+1;
        end
    end

    % Solve the roots of Legendre
    rxd=zeros(m-1,1); 
    for i=1:fix(m/2)
        if abs(rdy(i,1))<abs(rdy(i,2))
            rxd(i)=xdr(i,1);
        else
            rxd(i)=xdr(i,2);
        end
        for j=1:R
            [~,dyj,ddyj]=LegendreRecDer(m,rxd(i));
            rxd(i)=rxd(i)-dyj/ddyj;
        end
    end

    % Expand
    if mod(m,2)==0
        p=fix(m/2);
        for i=1:p-1
            rxd(p+i)=-rxd(p-i);
        end
    elseif mod(m,2)==1
        p=fix(m/2);
        for i=1:fix(m/2)
            rxd(p+i)=-rxd(p-i+1);
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


