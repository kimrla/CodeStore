function [uu, d, pp] =acrvints(ps1, ps2, ut)

% Get approximate intersection points of two set points
%
%  Inputs: 
%
%     ps1, ps2 - the two sets points to be used
% 
%     ut - parametric points of the two sets points
%
% Output: 
%
%     uu - the approximated intersection parametric points
%
%     pp - the approximated coordinates of the intersection 
%                            points of the two curves
% 
%     d - the minimum distances of the two curves

% Get the index of approximated nearest points according to distance
dm = DistanceMatrix(ps1',ps2'); 
[midm, II]=min(dm);
dmidm=mfdm(ut, midm);
n=length(ut); t=1; IJ=[]; 
for i=1:n-1
    if dmidm(i)<0 && dmidm(i+1)>0
        IJ(t,1)=II(i+1); 
        IJ(t,2)=i+1;
        t=t+1;
    end
end
if isempty(IJ)
    [~, i]=min(midm);
    IJ(t,1)=II(i); 
    IJ(t,2)=i;
end
uu=[ut(IJ(:,1))', ut(IJ(:,2))'];    

% Get the coordinates of approximated nearest points
[m,~]=size(IJ); d=zeros(m,1);
pp=zeros(3,2,m); 
for i=1:m
    d(i)=sqrt(sum((ps1(:,IJ(i,1))-ps2(:,IJ(i,2))).^2));
    pp(:,:,i)=[ps1(:,IJ(i,1)), ps2(:,IJ(i,2))];
end





