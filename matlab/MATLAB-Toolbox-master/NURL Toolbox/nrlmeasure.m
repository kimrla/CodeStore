function arclen = nrlmeasure(crv, k)

% Measure the length of a NURL curve
%
%  Input:
%    crv - a nurl curve
%    k - the number of intergration points in each interval
%         (default value is 10)
%
%  Output:
%   arclen - length of the curve
%

if nargin==1
    k=10;
end

[s, C]=GaussR(k,-1,1);
m=length(crv.intervals); 
arclen=0; s=s'; C=C';
for i=1:m-1
    ai=crv.intervals(i); bi=crv.intervals(i+1);
    tt=(bi-ai)*s/2+(bi+ai)/2;
    Ci=(bi-ai)*C/2;
    [~, dpntsi] = nrldeval (crv, tt);
    dr=zeros(1,k);
    for j=1:k
        dr(j)=norm(dpntsi(1:3,j));
    end
    arclen=arclen+sum(Ci.*dr);
end


