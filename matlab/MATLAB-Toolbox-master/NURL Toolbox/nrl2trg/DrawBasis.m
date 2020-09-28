% Draw a basis

function DrawBasis(v1,v2,Sk,S,T)

surf(S,T,Sk); hold on; 
shading interp; 
plotedge(S,T,Sk,'r');
[u, v]=meshgrid(v1, v2);
[u, v]=getst([u(:)'; v(:)'], 1);
plot(u(:), v(:), 'ro'); 
hold off;