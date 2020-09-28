function  [pntt, order]=sort4nodes(pnts)

% SORTNODES: Sort the 4 nodes for constructs a NURL bilinear surface.
% 
% INPUT: 
% 
%     pnts : a cell of pints, for example, 
%                 pnts={[0.0 0.0 0.3], [1.0 0.0 -0.3], [0.0 1.0 -0.3], [1.0 1.0 0.3]};
% 
% OUTPUT:  
% 
%     pntt :   the sorted 4 ponts 
%     order:  the new order of the 4 ponts 
%
% Description:
%     
%      This function preparing the 4 points for nrl4surf function.
%

pt=ppp(pnts);  
pt3=pnts; pt3(4)=[]; 
[~,ags]=piotg(pt3, pt);  
[~,J] = max(ags) ; 
agn=[3,1,2,4;   1,2,3,4;   2,1,3,4]; 
order=agn(J,:);
pntt=pnts; 
for j=1:4
    pntt{j}=pnts{agn(J,j)};
end



function [pt,pl]=ppp(pnts)

% PPP: Get preject of the 4th points on the plane of first 3 points
%             and the plane formed by the first 3 points
%
% INPUT: 
% 
%     pnts: a cell of pints, for example, 
%                 pnts={[0.0 0.0 0.3], [1.0 0.0 -0.3], [0.0 1.0 -0.3], [1.0 1.0 0.3]};
% 
% OUTPUT:  
% 
%       pt :  the coordinate of the prection of the 4th points 
%                on the plane of first 3 points
%       pl :  the parameters [A, B, C, D] of the plane formed by first 3
%                  points, the plane is Ax+By+Cz+D=0
%


% Get the plane of first 3 points
%   [A, B, C] is the normal vector of first 3 points
%   Ax+By+Cz+D=0 is equation of the plane
pln=zeros(3,4); pln(:,4)=1;
pln(1,1:3)=pnts{1};
pln(2,1:3)=pnts{2};
pln(3,1:3)=pnts{3};
A=det(pln(:,[2,3,4]));
B=-det(pln(:,[1,3,4]));
C=det(pln(:,[1,2,4]));
D=-det(pln(:,[1,2,3]));

x0=pnts{4}(1); y0=pnts{4}(2); z0=pnts{4}(3);

% The normal line of the plane passed the 4th point
t=-(A*x0+B*y0+C*z0+D)/(A^2+B^2+C^2);
x=x0+A*t;
y=y0+B*t;
z=z0+C*t;

pt=[x,y,z];
pl=[A,B,C,D];



function [yn,ags,vt]=piotg(pt3, pt)

% PIOTG: whether a point is inside or outside a triangular
% 
% INPUT: 
% 
%     pt3 : a cell of the coordinates of 3 points, for example, 
%                 pnts={[0.0 0.0 0.3], [1.0 0.0 -0.3], [0.0 1.0 -0.3], [1.0 1.0 0.3]};
%     pt:    the 4th point, for example, pt = [0.6, 0.5, -3.6];
% 
% OUTPUT:  
% 
%       yn :  yn = 1, if the point is outside a triangular the 3 points
%                    yn = 0, if the point is inside a triangular the 3 points
%       vt :    the 3 vectors from the pont pt to the 3 points pt3
%       ag :  the 3 angles between the 3 vectors of vt
%


vt=zeros(3,4);
for j=1:3
    vt(:,j)=pt3{j}-pt;     
end
vt(:,4)=vt(:,1);

aglvect = @(A, B) (180/pi)*acos(dot(A,B)/(norm(A)*norm(B)));
ags=zeros(3,1);
for j=1:3
    ags(j)=aglvect(vt(:,j+1), vt(:,j));
end
ag=sum(ags);

if ag<355
    yn=1;
else
    yn=0;
end




%% Demo
% pnts={[0.0 0.0 0.3], [1.0 0.0 -0.3], [0.0 1.0 -0.3], [1.0 1.0 0.3]};
% pnts={pnts{3}, pnts{1}, pnts{2}, pnts{4}};
% [pntt, order]=sort4nodes(pnts);
% srf1 = nrl4surf(pntt{1}, pntt{2}, pntt{3}, pntt{4});
% srf1 = nrlkntins(srf1,{0.5 0.5});
% figure; nrlaxisplot(srf1); 
% hold off;



