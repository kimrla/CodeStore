%% 中值滤波器
%x为原始信号，R为某点左右邻域大小.(窗口大小=2*R+1)
function [y] = myMedfilt(x, R)
y = x;
for (i = 1:length(x))
    if ((i+R)<= length(x) && (i-R)>= 1)
        BL = median(x((i-R):(i+R)));
    elseif ((i+R)<= length(x) && (i-R)< 1)
        BL = median(x(1:(i+R)));
    elseif ((i+R)> length(x) && (i-R)>= 1)
        BL = median(x((i-R):end));
    end
    y(i) = y(i)- BL;
end
end