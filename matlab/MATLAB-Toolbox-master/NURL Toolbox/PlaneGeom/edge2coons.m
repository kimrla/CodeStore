% Get edges for making coons surfaces from arbitrary edges
%        inputs: crvs - the given arbitrary edges
%       output: coons - edge numbers of each coons surface
%       examples: load testcoons crvs;
%                       coons=edge2coons(crvs);

function coons=edge2coons(crvs)

% Remove points
crvs=RePoints(crvs);

% Get the end points of the curves
n=numel(crvs);
pnts=zeros(3, 2*n);
for i=1:n
    pntx =  nrleval(crvs(i), [0, 1]);
    pnts(:, 2*i-1:2*i) = pntx;
end

% Get the duplicate points
dx=max(pnts(1,:))-min(pnts(1,:));
dy=max(pnts(2,:))-min(pnts(2,:));
dz=max(pnts(3,:))-min(pnts(3,:));
dist=max([dx, dy, dz]);
dm = DistanceMatrix(pnts',pnts'); 
id=cell(n,1); dp=zeros(n,1);
for i=1:2*n
    dd=find(dm(i,:)/dist<0.02);
    t=dd==i; dd(t)=[];
    id{i}=dd;
    dp(i)=length(id{i});
end

% Get the neighbouring curves of each curve
nbs=cell(n,1);
for i=1:n
    dd=ceil([id{2*i-1}, id{2*i}]/2);
    nbs{i}=dd;
end

% Get the quadrangle through neighbouring curves
coons=[];
for k=1:n
    ed=oppedge(k, nbs, pnts, dist, n);    
    coons=[coons; ed];
end
coons=sort(coons, 2);

% Remove duplicate elements
coons=RemDuplicate(coons);


function ed=oppedge(k, nbs, pnts, dist, n)

% Get the points of neighbouring curves
nbc=nbs{k}; 
m=length(nbc);
pntsk=pnts(:, 2*k-1:2*k);
nbpnts=zeros(3, 2*m);
for i=1:m
    p=nbc(i);
    nbpnts(:, 2*i-1:2*i)=pnts(:, 2*p-1:2*p);
end
dm = DistanceMatrix(pntsk', nbpnts'); 
d1 = dm(1, :)/dist<0.02;
d2 = dm(2, :)/dist<0.02;
dd=(d1 + d2)~=0;
nbpnts(:,dd)=[];


% Get the opposite dges
dm=DistanceMatrix(nbpnts', pnts'); 
dd = dm/dist<0.02;
sd=sum(dd);
de=[];
for i=1:n
    if sd(2*i-1) && sd(2*i)
        de=[de, i];
    end
end

% Get the quadrangle
p=length(de); ed=[];
for i=1:p
    dc=nbs{de(i)};
    dm=DistanceMatrix(dc', nbc'); 
    ds=[];
    for j=1:m
        dk=dm(:,j)==0;        
        ds=[ds, dc(dk)];
    end
    ds=[k, de(i), ds];
    ed=[ed; ds];
end


% Remove  points
function  crvs=RePoints(crvs)
n=numel(crvs);
pnts=[]; t=1; T=[];
for i=1:n
    if strcmp(crvs(i).form, 'POINTS')
        pnts=[pnts, crvs(i).coefs'];
        T(t)=i;
        t=t+1;
    end
end
crvs(T)=[];



%! Demo
% load testPlate sdata
% crvs=sdata.crvs;
% 
% % Get edges for coons surfaces
% coons=edge2coons(crvs);
% 
% % Get the coons surfaces
% [m,~]=size(coons);
% for i=1:m
%     nrbs={crvs(coons(i,1)), crvs(coons(i,2)), crvs(coons(i,3)), crvs(coons(i,4))};
%     tcrvs=trans4cnsrld(nrbs);
%     srfs(i)=nrbcoonsm(tcrvs{1}, tcrvs{2}, tcrvs{3}, tcrvs{4});
% end
% 
% % Plot
% figure; hold on;
% chs= nrbcrvplot(crvs, 50);
% axis equal;
% 
% % figure; hold on;
% for i=1:m
%     fh{i}=nrbsrfplot(srfs(i),[15,15]);
% end
% axis equal;



