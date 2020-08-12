function Njp_u = Njp(j, p , u, ui)
% 计算基函数Nj,p(u),ui为节点向量
 
if p == 0       % 0次B样条
    if (u >= ui(j)) && (u < ui(j+1))
        Njp_u = 1.0;
    else
        Njp_u = 0.0;
    end
else
    Length1 = ui(j+p) - ui(j);
    Length2 = ui(j+p+1) - ui(j+1);      % 支撑区间的长度
    if Length1 == 0.0       % 规定0/0 = 0
        Length1 = 1.0;
    end
    if Length2 == 0.0
        Length2 = 1.0;
    end
    Njp_u = (u - ui(j)) / Length1 * Njp(j, p-1, u, ui) ...
        + (ui(j+p+1) - u) / Length2 * Njp(j+1, p-1, u, ui);
end

