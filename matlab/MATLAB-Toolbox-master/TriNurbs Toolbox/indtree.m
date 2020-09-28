function totaltree=indtree(bndin,totaltree,startbnd,endbnd,parentbnd)

% Construct the tree data structure of bndin, with making startbnd as the
% root node  and endbnd the leaf node.

% Input:
%   bndin:M*2 matrix
%   startbnd:Root node of the tree
%   endbnd:Leaf node of the tree

% Output:
%   totaltree:The entire tree data structure


% l=length(bndin);
% if (~iscell(totaltree))
%     totaltree=[];
% end

% while (true)
% if (deter==0)
    [u,v]=find(bndin==startbnd);
   

    k=2./v;
    anotherstartbnd=bndin(u,k);    
    paid=anotherstartbnd==parentbnd;
    u(paid)=[];v(paid)=[];anotherstartbnd(paid)=[];
    numu=length(u);
    
    if (numu~=0)
        tem2=find(anotherstartbnd==endbnd);



        if (sum(tem2)~=0)
            totaltree=startbnd;
            return;

    
        else
            totaltree=[totaltree,startbnd];
                for i=1:numu
                    startbnd=anotherstartbnd(i);

                    totaltree=indtree(bndin,totaltree,startbnd,endbnd);
                    
                end
            
        end
    else
        
    
    
    
    
    
    
    
    
    
    
    
    
    
end




            
            
            


