function crvg=nrlglues(crvs)

% Glue several curves to form a single curve
% 
%  Inputs: 
%       crvs - the curves to be glued
%
%  Output: 
%       crvg - the glued curve
%       i - index, success (1) or failed (0)
%
% Example:  
%       n=5; 
%       crvs=nrlpolygon(n); 
%       crvg=nrlglues(crvs); 
%       nrlplot(crvg, 1000); 
% 

n=numel(crvs);
crvg=crvs(n);
crvs(n)=[]; n=n-1;
t=1; p=n+1; 
while n>0
    flag=zeros(n,1);
    for j=n:-1:1
        [crv, id]=nrlglue(crvg, crvs(j));
        if id
            crvg=crv;
            flag(j)=1;
        end
    end
    flag=flag==0;
    crvs=crvs(flag);
    n=numel(crvs);
    t=t+1;
    if t>p
        break;        
    end
end







