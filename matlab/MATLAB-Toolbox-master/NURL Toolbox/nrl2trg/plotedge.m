function plotedge(x,y,z,c)

[m,n]=size(x); 
r=x(:,1); s=y(:,1); t=z(:,1); plot3(r,s,t,c);
r=x(:,n); s=y(:,n); t=z(:,n); plot3(r,s,t,c);
r=x(1,:); s=y(1,:); t=z(1,:); plot3(r,s,t,c);
r=x(m,:); s=y(m,:); t=z(m,:); plot3(r,s,t,c);