function tsrf=tnrboptimize(tsrf, p2t2, pfix, dt, it, h0, k0, method)

% tnrboptimize: Optimize the triangles of a tri-nurbs surface by Newton's method
% 
% Calling Sequences:
% 
%       tsrf=tnrboptimize(tsrf, p2t2, pfix, dt, it, h0, k0, method)
% 
% 
%       tsrf - Triangular representation of a nurbs surface.
% 
%       p2t2 - The relations from points to triangles of tri-nurbs 
%                  surface 1. See also tnrbpts2tri.
%
%       pfix - Indexes of the fixed nodes.
%
%       dt  -  Time step.
%
%       h0  -  Mesh seed length.
%
%       k0  -  Force factor.
%
%       Method - The method used for interiation. The default method 1
%                is the forward Euler method. The method 2 is Velocity
%                Stormer Verlet method. 
%
% OUTPUT: 
%
%       tsrf - The tri-nurbs surface after Optimization.
%

if nargin==7
    method=1;
end

srf=tsrf.nurbs;
pop=true(tsrf.numbers(1),1);
pop(pfix)=false;
pop=sort(find(pop));
np=length(pop);
F1=zeros(np, 3);
if method==2
    V1=F1; F2=F1; V2=F1;        
end
dsrf=nrbderiv(srf);
X2=zeros(3,1);
for k=1:it
    for i=1:np
        X1=tsrf.points(pop(i), :);
        x1=tsrf.nodes(pop(i), :);
        F1(i,:)=k0*tnrbforce(tsrf, p2t2, pop(i), h0);        
        for j=1:3
            if method==1
                X2(j)=X1(j)+dt*F1(i,j);
            elseif method==2
                X2(j)=X1(j)+dt*V1(i,j)+F1(i,j)*dt^2/2;                
            end
        end
        [x2, X2]=nrbsrfreverse(srf, dsrf, x1, X2);
        tsrf.points(pop(i), :)=X2;
        tsrf.nodes(pop(i), :)=x2;
        if method==2
            F2(i,:)=k0*tnrbforce(tsrf, p2t2, pop(i), h0);
            V2(i,j)=V1(i,j)+(F1(i,j)+F2(i,j))*dt/2;
        end
    end
end

%% demo
% % The mesh seed length (h0) and force facter (k0)
% h0=1.8; k0=1;
% 
% % Create a plane surface
% af=0;
% crv1=nrbcirc(8, [0,0], af, pi-af);
% crv2=nrbcirc(8, [0,0], pi+af, 2*pi-af);
% crv2=nrbreverse(crv2);
% % srf=nrbruled(crv1, crv2);
% srf=nrbrevolve(crv1, [0,0,0], [1,0,0], pi);
% 
% % Transform a nurbs surface into triangular representation
% tsrf=nrb2tri(srf, h0);
% 
% figure; hold on; 
% tnrbplot(tsrf);
% axis equal; view(3);
% title('The surface before opimization.');
% 
% % Get the boundary of the tri-nurbs surface
% p2t2=tnrbpts2tri(tsrf);
% bp=triboundary(tsrf, p2t2);
% bp=bp{1};
% 
% % Optimize the triangles
% dt=0.1; it=10;
% pfix=sort(RemDuplicate(bp));
% tsrf=tnrboptimize(tsrf, p2t2, pfix, dt, it, h0, k0);
% 
% % Plot results
% figure; hold on; 
% tnrbplot(tsrf);
% plot3(tsrf.points(bp,1), tsrf.points(bp,2), tsrf.points(bp,3), 'r')
% axis equal; view(3);
% title('The surface after opimization.');




