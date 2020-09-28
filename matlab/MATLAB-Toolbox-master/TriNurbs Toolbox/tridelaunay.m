function tri=tridelaunay(n)

% tridelaunay: Delaunay triangulation of gridded nodes on a unit triangle.
%
% Calling Sequences:
% 
%       tri=tridelaunay(m)
% 
% INPUTS:
% 
%       N - The number of nodes on edges. All three edges have the same number.      
%
% OUTPUT:
% 
%     tri - Triangulation connectivity list, specified as an m-by-dim matrix, 
%            where m is the number of triangles or tetrahedra, and dim is the 
%            number of vertices per triangle. Each element in tri is a Vertex ID. 
%            Each row of tri contains the vertex IDs that define a triangle.
%
% Discritopn:
%  
%    See also: TrigNodeVect
% 

% Get the matrix
p=1; Num=nan(n); 
tri=zeros((n-1)^2,3);
for j=1:n
    for i=1:(n-j)+1
        Num(j,i)=p;
        p=p+1;
    end
end
p=1;
for j=1:n-1
    for i=1:(n-j)-1
        tri(p,1)=Num(j,i); tri(p,2)=Num(j,i+1); tri(p,3)=Num(j+1,i); 
        p=p+1;
        tri(p,1)=Num(j+1,i); tri(p,2)=Num(j,i+1); tri(p,3)=Num(j+1,i+1); 
        p=p+1;
    end
    if isempty(i)
        i=1;
    else
        i=i+1;
    end
    tri(p,1)=Num(j,i); tri(p,2)=Num(j,i+1); tri(p,3)=Num(j+1,i); 
    p=p+1;
end

%% demo
% % Create a triangle
% T=[0,0,0; 1,0,0; 1,1,1];
% 
% % Plot the triangle
% n=15;
% tri=tridelaunay(n);
% [u, v]=TrigNodeVect(n);
% Ps=tripoint(T(1,:), T(2,:), T(3,:), u, v);
% trisurf(tri,Ps(:,1),Ps(:,2),Ps(:,3));




