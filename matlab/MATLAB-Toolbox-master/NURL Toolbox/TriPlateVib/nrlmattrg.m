function [G, Gx, Gy, S, T, C, s, t]=nrlmattrg(tsrf)

% Get first derivatives of the nurl basis matrix in area coordinates of a triangle
% 
%  INPUT:
%
%    tsrf    -   a nurl triangle patch.
%
% OUTPUT:
% 
%    G - a matrix of the basis
% 
%    Gx, Gy - matrices of first derivatives with respect to area coordinates
% 
%    S, T - area coordinates transformed from natural coodinates tt
%
%    s, t - the natrual coordinates on a unit square domain
%
%  Description:
%
%     This function is similar as the nrltrgmat function. This function
%     limited in tensor product grids and also included integration nodes.
%     The nrltrgmat function can use scattered data and did not include
%     nodes.
%

% The total number of nodes
srf=tsrf.faces;
num1=srf.number(1); 
num2=srf.number(2); 
TN=num1*(num2-1)+1;

% Get first derivatives and weights of derivatives in Cartesian coordinates
[s, Cs]=GaussLobattoQ(srf.knots{1}, srf.order(1)+2);
[t, Ct]=GaussLobattoQ(srf.knots{2}, srf.order(2)+2);
[S, T, C]=getstc({s, t}, {Cs, Ct});
weights=squeeze(srf.coefs(4,:,:));
[pnti, jaci] = nrlgintvdeval (srf.order, weights, srf.knots, {s, t});
Gr=pnti'; Grx=jaci{1}'; Gry=jaci{2}'; 

% Sum the interpolation basis on the concentrated point
M=length(s); N=length(t);
G=zeros(M*N, TN); Gs=G; Gt=G;
pp=1:TN-1; qq=TN:TN-1+num1;
G(:,pp)=Gr(:,pp); G(:,TN)=sum(Gr(:,qq),2); 
Gs(:,pp)=Grx(:,pp); Gs(:,TN)=sum(Grx(:,qq),2); 
Gt(:,pp)=Gry(:,pp); Gt(:,TN)=sum(Gry(:,qq),2); 

% Transform basis to triangle (except the concentrated point)
Gx=0*G; Gy=Gx;
Ts=T(:); Ss=S(:);
pp=Ts~=1; qq=~pp;
c11=zeros(size(pp));
c21=c11; c22=c11;
c11(pp)=1./(1-Ts(pp));
c21(pp)=Ss(pp)./(1-Ts(pp)).^2;
c22(:)=1;
c11=repmat(c11, 1, TN);
c21=repmat(c21, 1, TN);
c22=repmat(c22, 1, TN);
Gx(pp,:)=c11(pp,:).*Gs(pp,:);
Gy(pp,:)=c21(pp,:).*Gs(pp,:)+c22(pp,:).*Gt(pp,:);

% Transform basis on the concentrated point to triangle
uu=M*(N-1)+1; vv=M*N;
gx=Gt(uu, :)-Gt(vv, :);
gy=Gt(uu, :);
rr=length(find(qq));
gx=repmat(gx, rr, 1);
gy=repmat(gy, rr, 1);
Gx(qq, :)=gx;
Gy(qq, :)=gy;


%% Demo - derivetives in area coordinates
% % Creat a triangular nurls patch
% crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
% af=pi/3; a=cos(af); b=sin(af);
% crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
% crv3 = nrlcirc(1, [0, 0], 0, af);
% srf = nrltrgcoons(crv1, crv2, crv3);
% 
% % Transform the surface into  a nurls triangle patch
% srf=nrlkntins(srf, {2, 1});
% srf = nrldegelev(srf, [2, 2]);
% tsrf=nrl2trg(srf); 
% srf=tsrf.faces;
% 
% figure; hold on;
% nrlplot(srf, [20, 20]);
% nrlaxisplot(srf);
% view(2); axis equal;
% 
% % Coefficients of blending function on triangles
% v1 = srf.knots{1}; v2 = srf.knots{2}; 
% num1=length(v1); num2=length(v2); 
% NB=num1+2*num2-3; 
% TN=num1*(num2-1)+1;
% 
% % Basis on an triangle
% [G, Gs, Gt, S, T, ~, s, t]=nrlmattrg(tsrf);
% m=length(s); n=length(t);
% 
% [s1, r1]=meshgrid(v2(1:end-1), v1(1:end));
% pts=[r1(:)', 0; s1(:)', 1];
% pnts=nrleval(srf, pts);
% p=pnts*G';
% dps=pnts*Gs';
% dpt=pnts*Gt';
% 
% % Test for gridded points
% figure; hold on;
% x=reshape(p(1,:), m, n);
% y=reshape(p(2,:), m, n);
% z=reshape(p(3,:), m, n);
% surf(x,y,z);
% quiver(p(1,:), p(2,:), dps(1,:), dps(2,:));
% quiver(p(1,:), p(2,:), dpt(1,:), dpt(2,:));
% colormap summer;      
% axis equal;
% axis([-0.1, 1.1, -0.1,0.95]);

%% Demo - interpolation in natrual coordinates
% % Order of the NURLS basis (order)
% order=3*[1, 1]; nodes=[16, 17];
% 
% % Create a sector
% crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
% af=pi/2; a=cos(af); b=sin(af);
% crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
% crv3 = nrlcirc(1, [0, 0], 0, af);
% srf = nrltrgcoons(crv1, crv2, crv3);
% srf = nrldegelev(srf, order-srf.order);
% iknts=nodes-srf.number;
% srf=nrlkntins(srf, {iknts(1), iknts(2)});
% 
% % Transform the nurls triangle patch for subsequent manipulations
% tsrf=nrl2trg(srf);
% figure; hold on;
% nrlctrlplot(tsrf.faces); view(2);
% nrlaxisplot(tsrf.faces);
% 
% % The total number of nodes
% srf=tsrf.faces;
% num1=srf.number(1); 
% num2=srf.number(2); 
% TN=num1*(num2-1)+1;
% 
% % Get first derivatives and weights of derivatives in Cartesian coordinates
% [G, Gx, Gy, S, T, C]=nrlmattrg(tsrf); 
% [M, N]=size(S); 
% 
% % Interpolation on triangles
% [s, t]=getst(srf.knots);
% s=s(:); t=t(:);
% s=s(1:TN); 
% t=t(1:TN); 
% z=funz(s,t);
% Z=G*z;
% Zx=Gx*z;
% Zy=Gy*z;
% 
% % Draw the basis
% [F, Fx, Fy]=funz(S, T);
% Z=reshape(Z, M, N); 
% Zx=reshape(Zx, M, N); 
% Zy=reshape(Zy, M, N); 
% figure; surf(S, T, Z); title('Z');
% figure; surf(S, T, Z-F); title('error of Z');
% 
% figure; surf(S, T, Zx); title('Zx');
% figure; surf(S, T, Zx-Fx); title('error of Zx');
% 
% figure; surf(S, T, Zy); title('Zy');
% figure; surf(S, T, Zy-Fy); title('error of Zy');

%% Demo - interpolation on Cartesian coordinates
% % Order of the NURLS basis (order)
% order=3*[1, 1]; nodes=[16, 17];
% 
% % Create a sector
% crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
% af=pi/2; a=cos(af); b=sin(af);
% crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
% crv3 = nrlcirc(1, [0, 0], 0, af);
% srf = nrltrgcoons(crv1, crv2, crv3);
% srf = nrldegelev(srf, order-srf.order);
% iknts=nodes-srf.number;
% srf=nrlkntins(srf, {iknts(1), iknts(2)});
% 
% % Transform the nurls triangle patch for subsequent manipulations
% tsrf=nrl2trg(srf);
% figure; hold on;
% nrlctrlplot(tsrf.faces); view(2);
% nrlaxisplot(tsrf.faces);
% 
% % The total number of nodes
% srf=tsrf.faces;
% num1=srf.number(1); 
% num2=srf.number(2); 
% TN=num1*(num2-1)+1;
% 
% % Get first derivatives and weights of derivatives in Cartesian coordinates
% [G, Gx, Gy, S, T, C, pnts]=nrlplanemattrg(tsrf);
% [M, N]=size(S); 
% 
% % Interpolation on triangles
% pnt=srf.coefs(1:3,:);
% x=pnt(1,1:TN)';
% y=pnt(2,1:TN)';
% z=funz(x,y);
% Z=G*z;
% Zx=Gx*z;
% Zy=Gy*z;
% 
% % Draw the results
% [F, Fx, Fy]=funz(pnts(1,:)', pnts(2,:)');
% X=reshape(pnts(1,:)', M, N);
% Y=reshape(pnts(2,:)', M, N);
% F=reshape(F, M, N);
% Fx=reshape(Fx, M, N);
% Fy=reshape(Fy, M, N);
% Z=reshape(Z, M, N); 
% Zx=reshape(Zx, M, N); 
% Zy=reshape(Zy, M, N); 
% figure; surf(X, Y, Z); title('Z');
% figure; surf(X, Y, Z-F); title('error of Z');
% 
% figure; surf(X, Y, Zx); title('Zx');
% figure; surf(X, Y, Zx-Fx); title('error of Zx');
% 
% figure; surf(X, Y, Zy); title('Zy');
% figure; surf(X, Y, Zy-Fy); title('error of Zy');


%% Demo - basis
% % Order of the NURLS basis (order)
% % The number of nodes (nodes), the basis to be draw (k) 
% k=2; order=3*[1, 1]; nodes=[6, 7];
% 
% % Create a sector
% crv1 = nrlline([0.0 0.0 0.0]',[1.0 0.0 0.0]');
% af=pi/3; a=cos(af); b=sin(af);
% crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
% crv3 = nrlcirc(1, [0, 0], 0, af);
% srf = nrltrgcoons(crv1, crv2, crv3);
% srf = nrldegelev(srf, order-srf.order);
% iknts=nodes-srf.number;
% srf=nrlkntins(srf, {iknts(1), iknts(2)});
% 
% % Transform the nurls triangle patch for subsequent manipulations
% tsrf=nrl2trg(srf);
% figure; hold on;
% nrlctrlplot(tsrf.faces); view(2);
% nrlaxisplot(tsrf.faces);
% 
% % The total number of nodes
% srf=tsrf.faces;
% num1=srf.number(1); 
% num2=srf.number(2); 
% TN=num1*(num2-1)+1;
% 
% % Get first derivatives and weights of derivatives in Cartesian coordinates
% [G, Gx, Gy, S, T, C]=nrlmattrg(tsrf); 
% [M, N]=size(S); 
% 
% % Draw the basis
% Zg=reshape(G(:, k), M, N); 
% Zgx=reshape(Gx(:, k), M, N); 
% Zgy=reshape(Gy(:, k), M, N); 
% figure; surf(S,T,Zg); hold; title('Zg');
% shading interp; 
% plotedge(S,T,Zg,'k');
% [y,x]=meshgrid(srf.knots{2}, srf.knots{1});
% x=x.*(1-y);
% plot3(x,y,0*y,'ro');
% 
% figure; surf(S,T,Zgx); hold; title('Zgx');
% shading interp; 
% plotedge(S,T,Zgx,'k');
% plot3(x,y,0*y,'ro');
% 
% figure; surf(S,T,Zgy); hold; title('Zgy');
% shading interp; 
% plotedge(S,T,Zgy,'k');
% plot3(x,y,0*y,'ro');
% 
% % Check the matrices
% pnts=srf.coefs(1:3,:);
% pnts=pnts(:,1:TN);
% ps=G*pnts';
% dpx=Gx*pnts';
% dpy=Gy*pnts';
% 
% figure; hold on;
% nrlctrlplot(srf); view(2);
% plot(ps(:,1), ps(:,2), 'ro');
% quiver(ps(:,1), ps(:,2), dpx(:,1), dpx(:,2));
% quiver(ps(:,1), ps(:,2), dpy(:,1), dpy(:,2));




