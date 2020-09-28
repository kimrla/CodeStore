%CONVEX_HULL
%INPUT:     Point Set:(n x 2)
%OUPUT:     HULL Point List: (x x 2)
%Samples:   convex_hull_demo( randn(200,2)*100)
 
function L=convex_hull_demo(P)
 
%Test Data
%data_size = data_size
%P = randi([-50,50],[data_size,2]);
[num,dimension] = size(P);
if dimension ~= 2
    error('dimension error')
end
 
 
P = sortrows(P,[1,2]);
 
%====Visual Lization
board_left = min(P(:,1));
board_right = max(P(:,1));
board_bottom = min(P(:,2));
board_up = max(P(:,2));
x_padding = (board_right- board_left)*0.1;
y_padding = (board_up- board_bottom)*0.1;
plot_range= [board_left - x_padding,board_right + x_padding,board_bottom-y_padding,board_up+y_padding];
 
clf;
scatter(P(:,1),P(:,2),'b.');
axis(plot_range);
hold on
%====VisualLization
 
%if there is only one or two point remain,return it
if num < 3
     L = P;
end
 
%STEP ONE: Upper Hull:
L_upper = P([1,2],:); %Take first two points
hull_handle = plot(L_upper(:,1),L_upper(:,2),'ob-');
for i = 3:num
    L_upper = [L_upper;P(i,:)]; %add the point into list
     
    while size(L_upper,1) >= 3
        l_size = size(L_upper,1);
        if side_check2(L_upper(l_size-2:l_size,:)) %Check if it is valid
            break;  %Quit if Valid
        else
            L_upper(l_size-1,:) = []; %Remove the inner point and continue if not
        end
        set(hull_handle,'XData',L_upper(:,1),'YData',L_upper(:,2));drawnow;
         
    end
    set(hull_handle,'XData',L_upper(:,1),'YData',L_upper(:,2));drawnow;
end
  
  
 %Visualization
 plot(L_upper(:,1),L_upper(:,2),'bo-');
 %Visualization
 
   
%STEP Two: Build the lower hull
L_lower = [P([num,num-1],:)]; % Add P(n) and P(n-1)
set(hull_handle,'XData',L_lower(:,1),'YData',L_lower(:,2));drawnow;
 
 
for i = num-2:-1:1
    L_lower = [L_lower;P(i,:)];
    while size(L_lower,1) >= 3
        l_size = size(L_lower,1);
       if side_check2(L_lower(l_size-2:l_size,:)) %Check if it is valid
            break;  %Quit if Valid
        else
            L_lower(l_size-1,:) = []; %Remove the inner point and continue if not
       end   
       set(hull_handle,'XData',L_lower(:,1),'YData',L_lower(:,2));drawnow;
    end
    set(hull_handle,'XData',L_lower(:,1),'YData',L_lower(:,2));drawnow;
end
 
L_lower([1,size(L_lower,1)],:) = [];
if isempty(L_lower)
    L = L_upper;
else
    L = [L_upper;L_lower(2:size(L_lower,1)-1,:)];
end
hold off;
return