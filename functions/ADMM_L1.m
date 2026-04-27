function [x, output] = ADMM_L1(A, b, pm)
%ADMM_L1  Solve the unconstrained L1 problem
%             min_x  0.5 * ||A x - b||_2^2  +  lambda * ||x||_1
%   via ADMM in scaled form.
%
%   Inputs
%     A        m-by-n sensing matrix
%     b        m-by-1 measurement vector
%     pm       struct of options:
%       pm.lambda    regularization parameter        (default 1e-5)
%       pm.delta     ADMM penalty parameter          (default 10*lambda)
%       pm.maxit     maximum iterations              (default 5*n)
%       pm.reltol    relative-change tolerance       (default 1e-6)
%       pm.x0        initial guess (n-by-1)          (default zeros)
%       pm.xg        ground truth, for diagnostics   (default x0)
%
%   Outputs
%     x        recovered signal
%     output   diagnostic struct:
%       output.relerr     relative change ||y - y_old|| / ||y||
%       output.obj        objective value at each iteration
%       output.res        ||Ax - b|| / ||b||
%       output.err        ||x - xg|| / ||xg||
%       output.time       cumulative wall-clock time
%       output.iter_count number of iterations performed
%
%   In this repository ADMM_L1 is used only as a warm start for the
%   truncated Huber method.

    [M, N] = size(A);
    start_time = tic;

    % --- parameters ------------------------------------------------------
    if isfield(pm, 'lambda'); lambda = pm.lambda; else; lambda = 1e-5;       end
    if isfield(pm, 'delta');  delta  = pm.delta;  else; delta  = 10*lambda;  end
    if isfield(pm, 'maxit');  maxit  = pm.maxit;  else; maxit  = 5*N;        end
    if isfield(pm, 'x0');     x0     = pm.x0;     else; x0     = zeros(N,1); end
    if isfield(pm, 'xg');     xg     = pm.xg;     else; xg     = x0;         end
    if isfield(pm, 'reltol'); reltol = pm.reltol; else; reltol = 1e-6;       end

    % --- factor (I + (1/delta) A A^T) once -------------------------------
    AAt = A * A';
    L   = chol(speye(M) + (1/delta) * AAt, 'lower');
    L   = sparse(L);
    U   = sparse(L');

    x   = zeros(N, 1);
    Atb = A' * b;
    y   = x0;
    u   = x;

    obj         = @(x) 0.5 * norm(A*x - b)^2 + lambda * norm(x, 1);
    output.pm   = pm;
    iter_count  = 0;

    for it = 1:maxit
        % x-update: soft thresholding
        x = shrink(y - u, lambda / delta);

        % y-update via Sherman-Morrison-Woodbury
        yold = y;
        rhs  = Atb + delta * (x + u);
        y    = rhs/delta - (A' * (U \ (L \ (A * rhs)))) / delta^2;

        % dual update
        u = u + (x - y);

        % diagnostics
        relerr   = norm(yold - y) / max([norm(yold), norm(y), eps]);
        residual = norm(A*x - b) / norm(b);

        output.relerr(it) = relerr;
        output.obj(it)    = obj(x);
        output.time(it)   = toc(start_time);
        output.res(it)    = residual;
        output.err(it)    = norm(x - xg) / norm(xg);
        iter_count        = iter_count + 1;

        if relerr < reltol && it > 2
            break;
        end
    end

    output.iter_count = iter_count;
end
