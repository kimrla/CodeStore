function [x, pnts1, pnts2, dt]=insertintersects(srf1, srf2, x, pnts1, pnts2, dt, crv1, crv2, pts1)

% Insert intersection points of two surfaces to intersects
% 
% Calling Sequences:
% 
%     [x, pnts1, pnts2, dt]=insertintersects(srf1, srf2, x, pnts1, pnts2, dt, crv1, crv2, pts1)
% 
% INPUTS:
%
%     srf1, srf2 - nurls surfaces 
%     x - parametric points of intersections
%     pnts1 - points of intersections on surface 1
%     pnts2 - points of intersections on surface 2
%     dt - distance of the two nearest points
%     crv1 - parametric curve 1
%     crv2 - parametric curve 2
%     pts1 - the parametric points to be used
% 
% OUTPUT:
% 
%     x - parametric points of intersections
%     pnts1 - points of intersections on surface 1
%     pnts2 - points of intersections on surface 2
%     dt - distance of the two nearest points
% 
%  See also:
%       srfsinterscts
% 

n=length(dt);
ds=zeros(n-1,1); 
for i=1:n-1
    ds(i)=norm(pts1(:,i+1)-pts1(:,i));
end
da=sum(ds)/(n-1);
p=ds>da; 
xi=[]; qi=zeros(n-1,1);
for i=1:n-1
    if p(i)
        ni=round(ds(i)/da);
        xn=linspace(crv1.knots(i), crv1.knots(i+1), ni+2); 
        xn([1, ni+2])=[];
        xi=[xi, xn];
        qi(i)=ni; 
    end
end
ptsi1=nrleval(crv1, xi); ptsi1=ptsi1(1:2,:);
ptsi2=nrleval(crv2, xi); ptsi2=ptsi2(1:2,:);
[pntsi1, jaci1]=nrldeval(srf1, ptsi1);
[pntsi2, jaci2]=nrldeval(srf2, ptsi2);
[xi, pntsi1, pntsi2, dti]=srfsintersct(srf1, srf2, 2, ptsi1, pntsi1, jaci1, ptsi2, pntsi2, jaci2);
nn=n+length(dti);
xn=zeros(4, nn);
pntsn1=zeros(3,nn);
pntsn2=pntsn1;
dtn=zeros(nn,1);
for i=1:n-1
    ni=sum(qi(1:i-1));
    s=i+ni;
    xn(:, s)=x(:, i);
    pntsn1(:, s)=pnts1(:, i); 
    pntsn2(:, s)=pnts2(:, i); 
    dtn(s)=dt(i); 
    if p(i)
        ii=ni+(1:qi(i)); 
        jj=s+1:s+qi(i); 
        xn(:, jj)=xi(:, ii); 
        pntsn1(:, jj)=pntsi1(:, ii); 
        pntsn2(:, jj)=pntsi2(:, ii); 
        dtn(jj)=dti(ii); 
    end
end
xn(:, end)=x(:, n);
pntsn1(:, end)=pnts1(:, n); 
pntsn2(:, end)=pnts2(:, n); 
dtn(end)=dt(n); 
x=xn;
pnts1=pntsn1; 
pnts2=pntsn2; 
dt=dtn; 










