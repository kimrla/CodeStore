function nkts=countknts(intervals, knots)

% Count the number of knots in each intervals
% 
%  Output:
% 
%     nkts - the number of knots in each intervals
%

m=length(intervals); 
nkts=zeros(1,m-1);
for i=1:m-1 
    a=intervals(i); b=intervals(i+1);
    p=knots>=a & knots<=b;
    k=length(find(p)); 
    nkts(i)=k;
end