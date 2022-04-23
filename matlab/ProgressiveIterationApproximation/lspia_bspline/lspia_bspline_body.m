% LSPIA per B-Spline
%
% parametri da impostare prima della chiamata a lspia_bspline_body:
% q_x, q_y              <- vettori punti di approssimazione
% q_tt                  <- parametri dei punti di approssimazione
% x, y                  <- vettore punti 
% tt                    <- parametri dei punti x, y
% tol1                  <- tolleranza sull'errore di approssimazione
% tol2                  <- tolleranza sull'avanzamento (diff. errori al
%                       passo k e al passo k-1)
% grafici               <- 1 per conclusivo, 3 per grafici su più passi
% knot_partition:
%                       1. per prendere i parametri dei punti
%                       2. per punti uniformi in tt(1)...tt(end)
%                       3. per punti  Whitney–Schoenberg
% g                     <- grado della B-Spline
% adjust_axis           <- 1 per svg

%B 样条的  LSPIA
%调用 lspia_bspline_body 之前要设置的 % 参数：
% q_x, q_y <- 逼近点向量
% q_tt <- 逼近点参数
% x, y <- 向量点
% tt <- 点 x, y 的参数
% tol1 <- 近似误差容限
% tol2 <- 进料公差（差异。误差在
% 步骤 k 和步骤 k-1)
% 图 <- 1 为结论性的，3 为多步图
% 结分区：
% 1.取点参数
% 2. 对于 tt (1) ... tt (end) 中的均匀点
% 3. 惠特尼 – 勋伯格积分
% g <- B 样条的度数
% 调整轴 <- 1 用于 svg


if (~exist('adjust_axis'))
    adjust_axis = 0;
end


p = [x;y];
a = min(tt);
b = max(tt);

% iterazioni da mostrare nel subplot (in ultima posizione è sempre mostrata
% l'ultima iterazione)
%在子图中显示的迭代百分比（始终显示最后一个位置最后一次迭代）

step_to_show = [0, 1, 3, 5, 7];

k = length(p)-g-1;
if (knot_partition == 1)
    t = [repmat(a,1,g+1), tt(3:end-2), repmat(b,1,g+1)];
elseif (knot_partition == 2)
    t = [repmat(a,1,g), linspace(a, b, k+2), repmat(b,1,g)];
elseif (knot_partition == 3)
    t = an_not_a_knot(g,tt);
end
bs = an_bspl(g, t, q_tt);%B样条基函数矩阵


btb = bs' * bs;

% Peso esatto
% autoval = eig(btb);
% mu = 2/(min(autoval)+max(autoval));

% Peso ottimizzato
% 准确重量
% autoval = eig (btb);
% mu = 2 / (min (autoval) + max (autoval));

% 优化权重

C = max(sum(btb, 2)); % prendo il massimo della somma di ogni riga
mu = 2/C; % definisco mu come da approssimazione
% 取每行总和的最大值
% 我将 mu 定义为近似值
err1 = tol1 + 1;
avanzamento = tol2 + 1;
k = 0;

if(grafici > 0)
    n_val = 5000;
    p_val = linspace(a, b, n_val);
    bs_val = an_bspl(g, t, p_val);
end
if(grafici == 3)
    n_grafico = 1;
    figure(2)
    subplot(3,2,n_grafico)
    plot((bs_val*p(1,:)')', (bs_val*p(2,:)')','r');
    hold on
    plot(x, y, 'go');
    
    if (adjust_axis)
        axis ij
        axis equal
    end
    plot(q_x, q_y, 'b.');
    title('k: 0')
end

% itero finché l'errore di approssimazione è superiore alla tolleranza e
% finché il passo di avanzamento è anch'esso superiore alla tolleranza
% 迭代直到逼近误差大于公差 e
% 只要进给步距也大于公差

while(err1 > tol1 && avanzamento > tol2)
    k = k+1;
    delta = [q_x - (bs*p(1,:)')'; q_y - (bs*p(2,:)')'];
    for(i=1:(n+1))
        adj(1, i) = mu * sum(bs(:, i)'.*delta(1, :));
        adj(2, i) = mu * sum(bs(:, i)'.*delta(2, :));
    end
    
    p = p + adj;
    
    for j=1:length(delta)
       norma(j) = norm(delta(:,j)) ;
    end
    tmp = max(norma);
    avanzamento = abs(err1 - tmp);
    err1 = tmp;
    
    if(k<10 || mod(k, 5)==0)
        fprintf("Errore al passo %1d: %5.5e\n", k, err1);
        %fprintf("Avanzamento al passo %1d: %5.5f\n", k, avanzamento);
    end
    
    if (grafici==3 && any(ismember(step_to_show, k)))
        n_grafico = n_grafico + 1;
        subplot(3,2,n_grafico)
        plot((bs_val*p(1,:)')', (bs_val*p(2,:)')','r');
        hold on
        plot(q_x, q_y, 'b.');
        
        title(strcat('k: ', num2str(k)))
        
        if (adjust_axis)
            axis ij
            axis equal
        end
    end
    
end
if(avanzamento <= tol2)
   fprintf("Interrotto per tolleranza raggiunta su avanzamento (%3.0e)\n", tol2); 
end


if (grafici > 0)
    figure(1)
    plot(q_x, q_y, 'b.', 'DisplayName', 'Punti di approssimazione');
    hold on
    plot((bs_val*p(1,:)')', (bs_val*p(2,:)')','r', 'DisplayName', 'Curva approssimata');
    plot(p(1,:), p(2,:), '*', 'DisplayName', 'Punti di controllo');
    legend
    
    if (adjust_axis)
        axis ij
        axis equal
    end
    
    plottitle = strcat('n-Q: ', num2str(n_q), ' - n: ', num2str(n), ' - k: ', num2str(k), ' - tol1: ', num2str(tol1), ' - tol2: ', num2str(tol2));
    title(plottitle)
    
    fprintf("Passo %1d. Errore: %5.5e.\n", k, err1);
end

if(grafici == 3)
    figure(2)
    n_grafico = n_grafico + 1;
    subplot(3,2,n_grafico)
    plot((bs_val*p(1,:)')', (bs_val*p(2,:)')','r');
    hold on
    plot(q_x, q_y, 'b.');
    %plot(x, y, 'go');
    title(strcat('k: ', num2str(k)))
    
    if (adjust_axis)
        axis ij
        axis equal
    end
end