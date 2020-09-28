function K=AssembleInP(K, Ke, Num)

% AssembleInP: Assembling element matrices for plane problem
% 
% Calling Sequences:
%
%     K=AssembleInP(K, Ke, Num)
%
% INPUTS:
%
%      K - global matrix before assemble
% 
%      Ke - element matrix
%
%      Num - node numbers
%
% OUTPUT:
%
%   K  -  global matrix after assemble
% 

N=length(Num);

% Assemble K
[~, f]=size(Ke);
switch f
    case 1
        for i=1:N
            p=Num(i); 
            for m=1:2
                K(2*(p-1)+m)=K(2*(p-1)+m)+Ke(2*(i-1)+m);
            end
        end
    otherwise
        for i=1:N
            p=Num(i); 
            for m=1:2
                for j=1:N
                    q=Num(j); 
                    for n=1:2
                        K(2*(q-1)+n, 2*(p-1)+m)=K(2*(q-1)+n, 2*(p-1)+m)+Ke(2*(j-1)+n, 2*(i-1)+m);
                    end
                end
            end
        end
end







