function tt=checktt(tt)

% Check whether tt is a cell array with row vectors
%
% Discription: the computation of basis take vectors
%    as row vectors so some subrotines need to check 
%    them. 

if iscell(tt)
    n=numel(tt);
    for i=1:n
        [m, ~]=size(tt{i});
        if m~=1
            tt{i}=tt{i}';
        end
    end
else
    [~, m]=size(tt);
    if m==1
        tt=tt';
    end
end






