function tsrf=SectorGeom(R, k, order, nodes)

% Create a sector plate NURL patch
af1=(k-1)*pi/2; af2=k*pi/2; 
a=R*cos(af1); b=R*sin(af1);
crv1 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
a=R*cos(af2); b=R*sin(af2);
crv2 = nrlline([0.0 0.0 0.0]',[a b 0.0]');
crv3 = nrlcirc(R, [0, 0], af1, af2);
srf = nrltrgcoons(crv1, crv2, crv3);
srf = nrldegelev(srf, order-srf.order);
iknts=nodes-srf.number;
srf=nrlkntins(srf, {iknts(1), iknts(2)});

% Transform the nurl triangle patch for subsequent manipulations
tsrf=nrl2trg(srf);