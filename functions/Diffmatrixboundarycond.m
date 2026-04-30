function [Lx] = Diffmatrixboundarycond(m, n, flag)
    if nargin < 3, flag = 1; end % 默认Neumann边界
    
    switch flag
        case 1 % Neumann边界条件（导数为零）
            L1 = spdiags([-ones(m,1), ones(m,1)], 0:1, m, m);
            L1(end, end) = 0;
            Lx = kron(speye(n), L1);
%             L2 = spdiags([-ones(n,1), ones(n,1)], 0:1, n, n);
%             L2(end, end) = 0;
%             Ly = kron(L2, speye(m));
            
        case 3 % 周期性边界
            rows = 1:m;
            cols_next = mod(rows, m) + 1;
            L1 = sparse(rows, rows, -ones(m,1), m, m) + ...
                 sparse(rows, cols_next, ones(m,1), m, m);
            Lx = kron(speye(n), L1);
%             rows = 1:n;
%             cols_next = mod(rows, n) + 1;
%             L2 = sparse(rows, rows, -ones(n,1), n, n) + ...
%                  sparse(rows, cols_next, ones(n,1), n, n);
%             Ly = kron(L2, speye(m));
            
        case 4 % Dirichlet边界（外部值为零）
            L1 = spdiags([-ones(m,1), ones(m,1)], 0:1, m, m);
            L1(end, end) = -1; % 最后一行为固定值
            Lx = kron(speye(n), L1);
%             L2 = spdiags([-ones(n,1), ones(n,1)], 0:1, n, n);
%             L2(end, end) = -1;
%             Ly = kron(L2, speye(m));
            
        case 5 % 反射边界
            L1 = spdiags([-ones(m,1), ones(m,1)], 0:1, m, m);
            L1(end, end-1) = 1; % 反射：u(end+1) = u(end-1)
            Lx = kron(speye(n), L1);
%             L2 = spdiags([-ones(n,1), ones(n,1)], 0:1, n, n);
%             L2(end, end-1) = 1;
%             Ly = kron(L2, speye(m));
    end
end