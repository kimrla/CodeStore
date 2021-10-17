function P1 = BLGTree(be,en,P,q,up)
if be+1==en
    P1 = P;
    return 
end
     P1 = P;
     C(1) = P1(en,2)-P1(be,2);
     C(2) = P1(be,1)-P1(en,1);
     C(3) = -C(2)*P1(be,2)-C(1)*P1(be,1);
     max = 0;
     ip = 0;
     q =q+1;
    for i = be:en
        distance = abs((C(1)*P1(i,1)+C(2)*P1(i,2)+C(3))/(sqrt(C(1)^2+C(2)^2)));
        if distance>max
            max = distance;
            ip =i;
        end
    end
    P1(ip,3) = q;
    P1(ip,4) = max;
    if P1(ip,4)>up
        P1(ip,4) = up;
    end
    up = min(P1(be,4),P1(ip,4));
    P1 = BLGTree(be,ip,P1,q,up);
    up = min(P1(ip,4),P1(en,4));
    P1 = BLGTree(ip,en,P1,q,up);
% else
%     max = 0;
%     for i = be:en
%         distance = sqrt((P(be,1)-P(i,1))^2+(P(be,2)-P(i,2))^2);
%         if distance>max
%             max = distance;
%             ip = i;
%         end
%     end
%     AA = max;
%     if AA>up
%         AA = up;
%     end
%     BB = ip;
%     BLGTree(be,ip,P(be,ip),AA);
%     BLGTree(ip,en,P(ip,en),AA);
end
