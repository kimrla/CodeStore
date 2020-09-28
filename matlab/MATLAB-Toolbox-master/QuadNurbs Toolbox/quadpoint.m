function [P, Jac, Hess]=quadpoint(quad, u, v)
%
% quadpoint: Get points on a quadrangle.
% 
% Calling Sequences:
% 
%     P=quadpoint(quad, u, v)
% 
%     [P, Jac, Hess]=quadpoint(quad, u, v)
% 
% INPUTS:
%
%      quad=[P1x, P1y, P1z
%                 P2x, P2y, P2z
%                 P3x, P3y, P3z
%                 P4x, P4y, P4z]
%           Four vertexes of a quadrangle.
%
%      u, v - Parametric coordinates of the quadrangle defined on [0, 1].
%
% OUTPUT:
% 
%     P - Points corresponds to the parametric coodinates.
%
%     Jac = [Pu, Pv]  -  Directions on the two parametric coodinates.
%
%     Hess = [Puu, Puv; Puv, Pvv]  -  Directions on the two parametric coodinates.
%

P1=quad(1,:); P2=quad(2,:); P3=quad(3,:); P4=quad(4,:); 
N1=(1-u).*(1-v); dN1du=-(1-v); dN1dv=-(1-u); dN1dudv=1; 
N2=u.*(1-v); dN2du=(1-v); dN2dv=-u; dN2dudv=-1; 
N3=u.*v; dN3du=v; dN3dv=u; dN3dudv=1; 
N4=(1-u).*v; dN4du=-v; dN4dv=(1-u); dN4dudv=-1; 
dim=length(P1);
% Use 4 vertices to construct the ruled surface
if dim==3
    x=P1(1)*N1+P2(1)*N2+P3(1)*N3+P4(1)*N4;
    y=P1(2)*N1+P2(2)*N2+P3(2)*N3+P4(2)*N4;
    z=P1(3)*N1+P2(3)*N2+P3(3)*N3+P4(3)*N4;
    P=[x(:), y(:), z(:)];
elseif dim==2
    x=P1(1)*N1+P2(1)*N2+P3(1)*N3+P4(1)*N4;
    y=P1(2)*N1+P2(2)*N2+P3(2)*N3+P4(2)*N4;
    P=[x(:), y(:)];
end
if nargout>1
    if dim==3
        xu=P1(1)*dN1du+P2(1)*dN2du+P3(1)*dN3du+P4(1)*dN4du;
        yu=P1(2)*dN1du+P2(2)*dN2du+P3(2)*dN3du+P4(2)*dN4du;
        zu=P1(3)*dN1du+P2(3)*dN2du+P3(3)*dN3du+P4(3)*dN4du;
        xv=P1(1)*dN1dv+P2(1)*dN2dv+P3(1)*dN3dv+P4(1)*dN4dv;
        yv=P1(2)*dN1dv+P2(2)*dN2dv+P3(2)*dN3dv+P4(2)*dN4dv;
        zv=P1(3)*dN1dv+P2(3)*dN2dv+P3(3)*dN3dv+P4(3)*dN4dv;
        Pu=[xu(:), yu(:), zu(:)];
        Pv=[xv(:), yv(:), zv(:)];
    elseif dim==2
        xu=P1(1)*dN1du+P2(1)*dN2du+P3(1)*dN3du+P4(1)*dN4du;
        yu=P1(2)*dN1du+P2(2)*dN2du+P3(2)*dN3du+P4(2)*dN4du;
        xv=P1(1)*dN1dv+P2(1)*dN2dv+P3(1)*dN3dv+P4(1)*dN4dv;
        yv=P1(2)*dN1dv+P2(2)*dN2dv+P3(2)*dN3dv+P4(2)*dN4dv;
        Pu=[xu(:), yu(:)];
        Pv=[xv(:), yv(:)];
    end
    Jac={Pu, Pv};
end
if nargout>2
    if dim==3
        xuv=P1(1)*dN1dudv+P2(1)*dN2dudv+P3(1)*dN3dudv+P4(1)*dN4dudv;
        yuv=P1(2)*dN1dudv+P2(2)*dN2dudv+P3(2)*dN3dudv+P4(2)*dN4dudv;
        zuv=P1(3)*dN1dudv+P2(3)*dN2dudv+P3(3)*dN3dudv+P4(3)*dN4dudv;
        Puv=[xuv(:), yuv(:), zuv(:)];
    elseif dim==2
        xuv=P1(1)*dN1dudv+P2(1)*dN2dudv+P3(1)*dN3dudv+P4(1)*dN4dudv;
        yuv=P1(2)*dN1dudv+P2(2)*dN2dudv+P3(2)*dN3dudv+P4(2)*dN4dudv;
        Puv=[xuv(:), yuv(:)];
    end
    Hess={0, Puv; Puv, 0};
end


%% Test
% Q=[0,0,0; 1,0,1; 1,1,0; 0,1,1];
% syms u v
% 
% % The bases
% N=[(1-u).*(1-v); u.*(1-v); u.*v; (1-u).*v];
% dNdu=diff(N, u);
% dNdv=diff(N, v);
% dNdudv=diff(dNdu, v);
% 
% % The points
% P=0; dPdu=0; dPdv=0; dPdudv=0;
% for i=1:4
%     P=P+Q(i,:)*N(i);
%     dPdu=dPdu+Q(i,:)*dNdu(i);
%     dPdv=dPdv+Q(i,:)*dNdv(i);
%     dPdudv=dPdudv+Q(i,:)*dNdudv(i);
% end
% 
% % Test of quadpoint
% m=5; n=6;
% s=linspace(0,1,m);
% t=linspace(0,1,n);
% [t, s]=meshgrid(t, s);
% [Pt, Jac, Hess]=quadpoint(Q, s, t);
% 
% figure; ezsurf(P(1), P(2), P(3), [0,1,0,1]);
% shading interp; hold on;
% plot3(Pt(:,1), Pt(:,2), Pt(:,3), 'r*');





