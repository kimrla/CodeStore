function trmparaplot(pcrv, nuv, nt)

% Plot the parametric region of a trimmed nurls surface
%
% Calling Sequence:
% 
%   trmparaplot(pcrv, nuv, nt)
% 
% INPUT:
% 
%    pcrv : parametric trimmed curve of the surface
% 
%   nuv : the number of sampling points on each direction [nu, nv]
%
%    nt  :  the number of parametric points of the trimmed curve
%   

t1=linspace(0, 1, nuv(1));
t2=linspace(0, 1, nuv(2));
tt={t1, t2};

[v1, u1]=meshgrid(tt{2}, tt{1});

tf = ishold;
plot(u1, v1); hold on;
plot(u1', v1');
nrlplot(pcrv, nt);

if tf
    hold on;
else
    hold off;
end





