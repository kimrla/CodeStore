function triplotunit(S, T, Z)

% triplotunit: Mesh plot of gridded nodes on a unit triangle.
%
% Calling Sequences:
% 
%       triplotunit(S, T)
% 
%       triplotunit(S, T, Z)
% 
% INPUTS:
% 
%       S, T, Z  -  The nodes on the triangle.
% 
% Discription:
%
%      See also: TrigNodeVect
%

n=(sqrt(8*length(S)+1)-1)/2;
tf = ishold;
p=1; Cs=nan(n); Ct=Cs;
if nargin==3
    Cz=Cs;
end
for j=1:n
    for i=1:(n-j)+1
        Cs(j,i)=S(p);
        Ct(j,i)=T(p);
        if nargin==3
            Cz(j,i)=Z(p);
        end
        p=p+1;
    end
end
if nargin==2
    plot(Cs,Ct); hold on;
    plot(Cs',Ct');
    for j=1:n
        s=zeros(n-j+1,1); t=s;
        for i=1:n-j+1
            k=n-j+1;
            s(i)=Cs(i, k-i+1);
            t(i)=Ct(i, k-i+1);        
        end
        plot(s,t);
    end
elseif nargin==3
    plot3(Cs,Ct,Cz); hold on;
    plot3(Cs',Ct',Cz');
    for j=1:n
        s=zeros(n-j+1,1); t=s; z=s;
        for i=1:n-j+1
            k=n-j+1;
            s(i)=Cs(i, k-i+1);
            t(i)=Ct(i, k-i+1);       
            z(i)=Cz(i, k-i+1); 
        end
        plot3(s,t,z);
    end
end
if tf
    hold on;
else
    hold off;
end

%% Demo
% % Generate nodes
% n=6;
% [S, T]=TrigNodeVect(n);
% 
% % Plot unit triangle
% figure; hold on;
% plot([0,1,0,0], [0,0,1,0]);
% plot(S, T, 'ro');
% axis equal;
% triplotunit(S,T,0*S);






