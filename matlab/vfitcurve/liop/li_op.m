function [N,L,A] = li_op(A,r,P,L)
for i = L:length(P)-1
    %%  判定该段线是否会与圆相交， 判断方法为该折线是否序号靠前的端点在圆的内部，另一端是否在圆的外部
    C(1) = P(i+1,2)-P(i,2);  %% 求出该段线方程系数
    C(2) = P(i,1)-P(i+1,1);
    C(3) = -C(2)*P(i,2)-C(1)*P(i,1);
    distance3 = abs((C(1)*A(1)+C(2)*A(2)+C(3))/(sqrt(C(1)^2+C(2)^2)));  %% 该折线到圆心的距离
% %     distance1 = sqrt((A(1)-P(i,1))^2+(A(2)-P(i,2))^2);
    distance2 = sqrt((A(1)-P(i+1,1))^2+(A(2)-P(i+1,2))^2);
%     distance4 = sqrt(distance1^2-distance3^2);
%     distance5 = sqrt(distance2^2-distance3^2);
%     distance6 = sqrt((P(i+1,1)-P(i,1))^2+(P(i+1,2)-P(i,2))^2);
%     &&((distance6 == distance4+distance5)
    %%  
    if (distance3<r)&&(distance2>r)
%         distance4 = sqrt((A(1)-P(i,1))^2+(A(2)-P(i,2))^2);
        if C(2) == 0
            x = -C(3)/C(1);
            y1 = sqrt(r^2-(x-A(1))^2)+A(2);
            y2 = -sqrt(r^2-(x-A(1))^2)+A(2);
            distance7 = sqrt((P(i+1,1)-x)^2+(P(i+1,2)-y1)^2);
            distance8 = sqrt((P(i+1,1)-x)^2+(P(i+1,2)-y2)^2);
            if distance7>distance8
                y = y2;
                N = [real((A(1)+x)/2),real((A(2)+y)/2)];
                A = [x,y];
                L = i;
            else
                y = y1;
                N = [real((A(1)+x)/2),real((A(2)+y)/2)];
                A = [x,y];
                L = i;
            end
        else
        a = sqrt(1+(C(1)/C(2))^2);
        a2b = (C(1)*C(3))/(C(2)^2)-A(1)+(C(1)*A(2))/C(2);
        b = a2b/a;
        x1 = real((sqrt(r^2-A(1)^2-(C(3)/C(2))^2-(2*A(2)*C(3))/C(2)-A(2)^2+b^2)-b)/a);
        y1 = real(-(C(1)*x1+C(3))/C(2));
        x2 = real((sqrt(r^2-A(1)^2-(C(3)/C(2))^2-(2*A(2)*C(3))/C(2)-A(2)^2+b^2)+b)/(-a));
        y2 = real(-(C(1)*x2+C(3))/C(2));
        distance7 = sqrt((P(i+1,1)-x1)^2+(P(i+1,2)-y1)^2);
        distance8 = sqrt((P(i+1,1)-x2)^2+(P(i+1,2)-y2)^2);
        if distance7>distance8
            x = x2;
            y = y2;
        else
            x = x1;
            y = y1;
        end
        N = [real((A(1)+x)/2),real((A(2)+y)/2)];
        L = i;
        A = [x,y];
        break
        end
    end
end




% r = 1.1;
% P(1,:) = [2,0];
% P(2,:) = [1,2];
% A = [1,1];
% C(1) = P(2,2)-P(1,2);
% C(2) = P(1,1)-P(2,1);
% C(3) = -C(2)*P(1,2)-C(1)*P(1,1);
% a = sqrt(1+(C(1)/C(2))^2);
% a2b = (C(1)*C(3))/(C(2)^2)-A(1)+(C(1)*A(2))/C(2);
% b = a2b/a;
% x(:) = (sqrt(r^2-A(1)^2-(C(3)/C(2))^2-(2*A(2)*C(3))/C(2)-A(2)^2+b^2)-b)/a;
% y(:) = -(C(1)*x+C(3))/C(2);