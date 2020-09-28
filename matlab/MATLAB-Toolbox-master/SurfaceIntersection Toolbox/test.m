clear; clc;

% Create Surface #1: the double cone
XYZ=[];
XYZ(1,:)=[0 0 0];
for z = 1:100
  n=z+5;
  [x,y] = pol2cart((1:n)'*2*pi/n,z);
  XYZ = [XYZ; [x y z*x./x]];
end
DT = delaunayTriangulation(XYZ(:,1),XYZ(:,2));
Surface1.faces = [DT.ConnectivityList; DT.ConnectivityList+size(XYZ,1)];
Surface1.vertices = [XYZ(:,3),XYZ(:,2),XYZ(:,1); -XYZ(:,3),XYZ(:,2),XYZ(:,1)]/100;
S=Surface1; trisurf(S.faces, S.vertices(:,1),S.vertices(:,2),S.vertices(:,3),'FaceAlpha', 0.5, 'FaceColor', 'r');
% Create Surface #2: 11 equally spaced parallel planes
Surface2=[];
for i=0:10
  z = -0.95 + i/5;
  Surface2.vertices(3*i+(1:3),:) = [2, 0, z; -1, 1.7, z; -1, -1.7, z];
  Surface2.faces(i+1,:) = 3*i+(1:3);
end
[~, Surf12] = SurfaceIntersection(Surface1, Surface2);
Surf12.vertices(:,3) = -1.5; % project the contour lines on a single plane
% plot the results
figure(1); clf; hold on
S=Surface1; trisurf(S.faces, S.vertices(:,1),S.vertices(:,2),S.vertices(:,3),'FaceAlpha', 0.5, 'FaceColor', 'r')
S=Surface2; trisurf(S.faces, S.vertices(:,1),S.vertices(:,2),S.vertices(:,3),'FaceAlpha', 0.5, 'FaceColor', 'b')
S=Surf12;   line(...
  [S.vertices(S.edges(:,1),1), S.vertices(S.edges(:,2),1)]',...
  [S.vertices(S.edges(:,1),2), S.vertices(S.edges(:,2),2)]',...
  [S.vertices(S.edges(:,1),3), S.vertices(S.edges(:,2),3)]',...
  'Color', 'r');
axis equal
view(-190, 15)




