%% ========================================================================
function PlotTwoTriangles(points, polygon, color)
% Plotting function used for debugging
nx = size(points,1)-6;
d = (max(points,[],1)-min(points,[],1))/200;
figure(2)
clf
hold on
line( points(nx+(1:2),1), points(nx+(1:2),2), points(nx+(1:2),3), 'Color', 'g');
line( points(nx+(2:3),1), points(nx+(2:3),2), points(nx+(2:3),3), 'Color', 'g');
line( points(nx+[1,3],1), points(nx+[1,3],2), points(nx+[1,3],3), 'Color', 'g');
line( points(nx+(4:5),1), points(nx+(4:5),2), points(nx+(4:5),3), 'Color', 'b');
line( points(nx+(5:6),1), points(nx+(5:6),2), points(nx+(5:6),3), 'Color', 'b');
line( points(nx+[4,6],1), points(nx+[4,6],2), points(nx+[4,6],3), 'Color', 'b');
plot3( points(:,1), points(:,2), points(:,3), 'm.');
if (length(polygon)>2)
  idx = polygon([1:end, 1]);
  plot3( points(idx,1), points(idx,2),points(idx,3), 'Color', color, 'LineWidth', 1);
end
for i = 1:nx+6
  text(points(i,1)+d(1), points(i,2)+d(2), points(i,3), num2str(i))
end

end % function