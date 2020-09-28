function nurbs=nrl2nrb(nurls)

% Transform  a nurl curve to a nurbs curve

p=nurls.order;
m=length(nurls.intervals);
coefs=zeros(4, (m-1)*p+1);
knots=zeros(1, m*p+2); 
knots(end-p:end)=1;
ni=2;
for i=1:m-1
    a=nurls.intervals(i); b=nurls.intervals(i+1);     
    Ui=[a*ones(1, p+1), b*ones(1, p+1)];
    ui=linspace(a, b, ni+p-1);
    si = findspan (ni, p, ui, Ui); 
    Ni = numbasisfun (si, ui, p, Ui)+1;
    Bi = basisfun (si, ui, p, Ui);
    Ci=zeros(ni+p-1);
    for j=1:ni+p-1
        Ci(j, Ni(j,:))=Bi(j, :);
    end
    [pwi, wi]=nrleval(nurls, ui);
    coefsi=(Ci\[pwi; wi]')';
    coefs(:, (i-1)*p+1 : i*p+1) = coefsi(:, :);
    knots((i-1)*p+2 : (i-1)*p+p+1)=a;
end

nurbs=nrbmak(coefs, knots);


%% Demo
% p=2; dt=1; 
% x=[0.5 1.5 4.5 3.0 7.5 6.0 8.5];
% y=[3.0 5.5 5.5 1.5 1.5 4.0 4.5];
% crv=nrlspline(p, x, y); 
% nrlctrlplot(crv);
% 
% crvb=nrl2nrb(crv);
% figure;
% nrbctrlplot(crvb);





