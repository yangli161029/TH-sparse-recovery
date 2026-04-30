function [unew,w,result] = algvetgradiffbdy(u, lambda, thr, maxit, flag) % 新增flag参数
    [m,n] = size(u);
    g = u;
    RE = 1e-8;
    i_list = [];
    
    for outloop = 1:40
        i = 0;
        for k = 1:maxit
            % 根据flag选择ux计算方式
            if flag == 1 % Neumann边界条件（导数为零）
                ux = zeros(size(u));
                ux(1:end-1, :, :) = u(2:end, :, :) - u(1:end-1, :, :);
            elseif flag == 3 % 周期性边界
                ux = u([2:end 1], :, :) - u;
            elseif flag == 4 % Dirichlet边界（外部值为零）
                ux = zeros(size(u));
                ux(1:end-1, :, :) = u(2:end, :, :) - u(1:end-1, :, :);
                ux(end, :, :) = -u(end, :, :); % 假设外部值为0
            elseif flag == 5 % 反射边界
                ux = zeros(size(u));
                ux(1:end-1, :, :) = u(2:end, :, :) - u(1:end-1, :, :);
                ux(end, :, :) = u(end-1, :, :) - u(end, :, :);
            end
            
            % 权重计算与原代码一致
            w = ones(size(ux));
            w(abs(ux) > thr) = 0;
            
            % 生成差分矩阵（需同步修改SparseDiffMatrix）
            [L] = Diffmatrixboundarycond(m, n, flag); % 传递flag参数
            D = spdiags(w(:), 0, m*n, m*n);
            L = L' * D * L;
            
            % 后续求解与原代码一致
            unew = vecsubprb(g, lambda, L);
            relerror = norm(unew - u) / norm(unew);
            i = i + 1;
            if relerror < RE
                break;
            end
            u = unew;
        end
        i_list = [i_list, i];
        thr = thr / 1.225;
        if thr < 1e-8
            break;
        end
    end
    result.i_list = i_list;
end