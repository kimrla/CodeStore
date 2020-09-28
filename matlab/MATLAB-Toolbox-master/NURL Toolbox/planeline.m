function varargout = planeline(varargin)

% planeline: Create lines and relations of plane lines through analytical expressions.
%       '2points' : create a line by two points
%       'pointtang' : create a line by a point and tangent of the line
%       'pointnormal' : create by a point and normal of the line
%       'intersection' : intersection and angle of two lines
%       'distance': distance of point to line
%       'center3ps':  get the center of the triangle formed by three points
%  输出：二维平面直线ABC:[a,b,c]:ax+by+c=0
% Case '2points' : create a line by two points
%           inputs : p1 - point 1, p2 - point 2
%       example : pnt1=[1,0];  pnt2=[2,2]; 
%                       ABC= planeline('2points',pnt1,pnt2);
% Case 'pointtang' : create a line by a point and tangent of the line
%           inputs : p1 - point 1, tang - tangent vector of the line
%       example : pnt1=[1,0];  tang=[1,-2]; 
%                       ABC= planeline('pointtang',pnt1,tang);
% Case 'pointnormal' : create a line by a point and normal of the line
%           inputs : p1 - point 1, norm - normal vector of the line
%       example : pnt1=[1,0];  norm=[1,2]; 
%                       ABC= planeline('pointnormal', pnt1, norm);
% Case 'intersection' : intersection and angle of two lines
%           inputs : ABC1 - line 1, ABC2 - line 2
%       example : ABC1=[-1,-2,3.5];  ABC2=[2,1,-3.5]; 
%                       [inter,angle]=planeline('intersection',ABC1,ABC2 );
% Case 'distance': distance of point to line
%            inputs : pnt - point, ABC - line
%       example : pnt=[0.5,0];  ABC=[-0.5,-1,0.875]; 
%                       d=planeline('distance', pnt, ABC );
% Case 'center3ps':  get the center of the triangle formed by three points
%           inputs : p1 - point 1, p2 - point 2, p3 - point 3
%       example : pnt1=[0.5,0]; pnt2=[1,1]; pnt3=[0,1]; 
%                       center= planeline('center3ps',pnt1,pnt2,pnt3);

switch varargin{1}
    case '2points'
        pnt1=varargin{2};
        pnt2=varargin{3};
        ABC(1)=det([pnt1(2), 1; pnt2(2), 1;]);
        ABC(2)=-det([pnt1(1), 1; pnt2(1), 1;]);
        ABC(3)=det([pnt1(1), pnt1(2); pnt2(1), pnt2(2);]);
        varargout{1}=ABC;
    case 'pointtang'
        pnt1=varargin{2};
        tang=varargin{3};
        pnt2=pnt1+tang;
        ABC= planeline('2points',pnt1,pnt2 );
        varargout{1}=ABC;
    case 'pointnormal'
        pnt1=varargin{2};
        norm=varargin{3};
        tang(1)=-norm(2);
        tang(2)=norm(1);
        pnt2=pnt1+tang;
        ABC= planeline('2points',pnt1,pnt2 );
        varargout{1}=ABC;
    case 'intersection'
        ABC1=varargin{2};
        ABC2=varargin{3};
        A1=ABC1(1); B1=ABC1(2); C1=ABC1(3); 
        A2=ABC2(1); B2=ABC2(2); C2=ABC2(3); 
        pp=A1*B2-A2*B1;
		mm=max(abs([A1,A2,B1,B2]))*1e-6;
        if abs(pp)>mm
            inter(1)=(B1*C2-B2*C1)/pp;
            inter(2)=(C1*A2-C2*A1)/pp;
            qq=sqrt(A1^2+B1^2)*sqrt(A2^2+B2^2);
            angle=acos((A1*A2+B1*B2)/qq);
        else
            inter=nan;
            angle=0;
        end
        varargout{1}=inter;
        varargout{2}=angle;
    case 'distance'
        pnt=varargin{2};
        ABC=varargin{3};
        A=ABC(1); B=ABC(2); C=ABC(3); 
        d=abs(A*pnt(1)+B*pnt(2)+C)/sqrt(A^2+B^2);
        varargout{1}=d;
    case 'center3ps'
        pnt1=varargin{2};
        pnt2=varargin{3};
        pnt3=varargin{4};
        
        % Middle points
        pm1=(pnt1+pnt2)/2; 
        pm2=(pnt2+pnt3)/2; 

        % Tangent vector 
        tng1=pnt2-pnt1;
        tng2=pnt3-pnt2;

        % Normal line equations
        nlin1= planeline('pointnormal',pm1,tng1 );
        nlin2= planeline('pointnormal',pm2,tng2 );

        % Get the intersection ponit and angle of two lines
        center=planeline('intersection', nlin1, nlin2 );
        
        varargout{1}=center;
end






