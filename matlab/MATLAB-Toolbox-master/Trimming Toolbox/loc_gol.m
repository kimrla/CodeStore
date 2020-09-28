function [nuv,mul]=loc_gol(bzr,uv,mode)

% Transform the par-coords in global domain to that in local domian, or
% vice versa. 

% Input:
%   bzr: Bezier structure.
%   uv: Par-coords in global/local domain, which can be more than one:
%        for curve:[u1,u2,...] or for surface:{[u1,u2,...],[v1,v2,...]}.
%   mode: 'local': From global u to local t; 'global': From local t to global u.
% Output:
%   nuv: Par-coords in local/global domain, which corresponds to the input uv. nuv is a vector(for curve) or a cell array(for surface).
%   mul: Multiple between the first order tangent vector of the former
%       domian to the latter domain. 'local': dp/du=mul*dp/dt; 'global':
%       dp/dt=1/mul*dp/du. mul=b-a or 1/(b-a).

knots=bzr.knots;
% num=length(uv);
nuv=[];
if iscell(knots)% Surface
    a(1)=min(knots{1});a(2)=min(knots{2});
    b(1)=max(knots{1});b(2)=max(knots{2});
else% curve
    a=min(knots);b=max(knots);
    uv=mat2cell(uv);
end
if strcmp(mode,'local')
    trans=@(u,a,b)((1-u).*a+u.*b);% Transform the global domain u=[0,1] to the local domain t=[a,b]
    mul=1./(b-a);
elseif strcmp(mode,'global')    
    trans=@(t,a,b)((t-a)./(b-a));% Transform the the local domain t=[a,b] to the global domain u=[0,1]
    mul=b-a;
end

if iscell(knots)
%     num=length(knots{1});
%     for i=1:num
%         nuv{1}=trans(uv{1},a(1),b(1));
%     end
%     num=length(knots{2});
%     for i=1:num
%         nuv{2}=trans(uv{2},a(2),b(2));
%     end
    nuv=cellfun(trans,uv,num2cell(a),num2cell(b),'UniformOutput',false);

else
%     num=length(knots);
%     for i=1:num
        nuv=trans(uv,a,b);
%     end
end
end


