function srf = nrlsphere(R, center, s1, s2, t1, t2)

% Create a nurl sphere
%
%  Input:
%    R - Radius
%    center - center
%    s1, s2 - start and end angle of x-y plane
%    t1, t2 - start and end angle of z-x plate
% 
%  Output:
%    srf - a nurl surface
%

% Creat x-y and z-y curves
crv1 = nrlellip(1, 1, [0, 0], s1, s2);
crv2 = nrlellip(1, 1, [0, 0], t1, t2);
crv2.coefs(3,:)=crv2.coefs(1,:);
crv2.coefs(1,:)=crv2.coefs(2,:);
crv2.coefs(2,:)=0;

% Make a nurl surface from the two curves
m=crv1.number; n=crv2.number;
order=[crv1.order, crv2.order];
knots={crv1.knots, crv2.knots};
intervals={crv1.intervals, crv2.intervals};
coefs=zeros(4, m, n);
for i=1:m
    for j=1:n
        coefs(1,i,j)=crv1.coefs(1,i)*crv2.coefs(1,j);
        coefs(2,i,j)=crv1.coefs(2,i)*crv2.coefs(1,j);
        coefs(3,i,j)=crv2.coefs(3,j);
        coefs(4,i,j)=crv1.coefs(4,i)*crv2.coefs(4,j);
    end
end
coefs(1,:,:)=R*coefs(1,:,:);
coefs(2,:,:)=R*coefs(2,:,:);
coefs(3,:,:)=R*coefs(3,:,:);
srf=nrlmake(coefs, knots, intervals, order);

trans = vectrans(center);
srf = nrltform(srf, trans);

%% demo
% R=1; N=6;
% s1=0; s2=pi; t1=0; t2=pi/2;
% center=[0, 0, 0];
% srf = nrlsphere(R, center, s1, s2, t1, t2);
% nrlplot(srf, [100, 100], 'ctrl');
% axis equal;



