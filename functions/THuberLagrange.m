function [x, result] = THuberLagrange(A, b, Param)
%THUBERLAGRANGE  Truncated Huber penalty for sparse signal recovery.
%
%   Implements the BCD algorithm (Algorithm 3.1 / 3.2) from
%       L. Yang, S. Morigi, M. K. Ng, Y.-W. Wen,
%       "Truncated Huber Penalty for Sparse Signal Recovery with
%        Convergence Analysis," SIAM J. Sci. Comput., 48 (2026),
%        pp. A929-A957.   DOI: 10.1137/25M1748184
%
%   Inputs
%     A          m-by-n sensing matrix (m << n)
%     b          m-by-1 measurement vector
%     Param      struct of options:
%       Param.x0        initial guess (n-by-1), required
%       Param.num       continuation: index into sorted |x| used to
%                       initialize mu (i.e. mu_0 = sort(|x0|,'descend')(num))
%       Param.DataRegP  Lagrangian-style fidelity weight (= alpha)
%       Param.xg        ground truth (only stored in result; unused inside)
%       Param.iter      max inner iterations per mu-epoch (default 100)
%       Param.RE        relative-change stopping tolerance (default 1e-8)
%
%   Outputs
%     x          recovered n-by-1 signal
%     result     struct with iteration count and a copy of inputs
%
%   The outer loop performs the mu-continuation strategy of Section 3.5;
%   the inner loop is the closed-form BCD update for fixed mu.

    iter     = 100;
    result   = Param;
    [M, N]   = size(A);
    dataregp = Param.DataRegP;
    RE       = 1e-8;

    if isfield(Param, 'RE');   RE   = Param.RE;   end
    if isfield(Param, 'iter'); iter = Param.iter; end

    %% --- Initialization -------------------------------------------------
    x      = Param.x0;
    tmp    = sort(abs(x), 'descend');
    num    = Param.num;
    mu     = tmp(num);

    idx        = abs(x) > mu;       % indicator of "large" entries
    sidx       = sum(idx);
    iter_count = 0;

    %% --- Outer loop: mu-continuation -----------------------------------
    for muiter = 1:40

        % Inner BCD loop with mu fixed
        for k = 1:iter
            A1    = A(:,  idx);
            A2    = A(:, ~idx);
            A2A2T = A2 * A2';
            HH    = inv(eye(M) / dataregp + (mu^2 / 2) * A2A2T);
            B     = A1' * (HH * A1);

            x1    = B \ (A1' * (HH * b));
            y     = HH * (A1 * x1 - b);
            x2    = -(mu^2 / 2) * (A2' * y);

            xnew         = x;
            xnew( idx)   = x1;
            xnew(~idx)   = x2;

            idx        = abs(xnew) > mu;
            sidx       = sum(idx);
            iter_count = iter_count + 1;

            if norm(xnew - x) / norm(xnew) < RE
                break;
            end
            x = xnew;
        end

        % --- Update mu (continuation rule, eq. (3.17)) -----------------
        muold1 = mu;
        m_I    = M - sidx;
        patl   = x(~idx);
        sum_p  = sum(patl(:).^2);

        % rho = 2.5 by default (Gaussian); use 1.7 for the oversampled DCT
        % SNR sweep, see comments in the original implementation.
        mu1   = muold1 / 2.5;
        mu2   = sqrt(sum_p / m_I);
        mu    = max([mu1, mu2]);

        idx   = abs(x) > mu;

        if mu < 1e-8, break; end
    end

    result.x          = x;
    result.iter_count = iter_count;
end
