function varargout = planepls(varargin)

% Get relations between planes and plane lines
% 
% Case 'pointnormal' : create a plane by point and normal of the plane
%           inputs : P - point, N - normal vector
%       example : P=[1, 0, 0];  N=[0, 0, 1]; 
%                       ABCD= planepls('pointnormal', P, N);
% 
% Case 'intersection' : intersection and angle of two lines
%           inputs : P1 - a point on line 1, T1 - direction of line 1
%                       P2 - a point on line 2, T2 - direction of line 2
%           outputs:
%                       d - the nearest distance of two curves
%                       P - coordinates on line 1
%                       Q - coordinates on line 2
%       example : P1=[0, 0, 0];  T1=[1, 0, 0]; 
%                       P2=[0, 0, 0];  T2=[0, 1, 0]; 
%                       [d, P, Q]=planepls('neardist', P1, T1, P2, T2);

switch varargin{1}
    case 'pointnormal'
        P0=varargin{2};
        N0=varargin{3};
        D=-(N0(1)*P0(1)+N0(2)*P0(2)+N0(3)*P0(3));
        varargout{1}=[N0(1), N0(2), N0(3), D];
    case 'neardist'
        P0=varargin{2};
        T0=varargin{3};
        P2=varargin{4};
        T2=varargin{5};
        X=P0-P2;
        A=zeros(2); B=zeros(2,1);
        for i=1:3
            A(1,1)=A(1,1)+T0(i)^2; 
            A(1,2)=A(1,2)-T0(i)*T2(i); 
            A(2,1)=A(2,1)-T0(i)*T2(i); 
            A(2,2)=A(2,2)+T2(i)^2; 
            B(1)=B(1)+T0(i)*X(i);
            B(2)=B(2)-T2(i)*X(i);
        end
        S=-A\B;
        P=P0+T0*S(1); Q=P2+T2*S(2);
        d=norm(P-Q);
        varargout{1}=d;
        varargout{2}=P;
        varargout{3}=Q;
end





