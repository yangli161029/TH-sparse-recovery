function A = data_generator_A(pm)
%DATA_GENERATOR_A  Build a sensing matrix used in the experiments.
%
%   A = DATA_GENERATOR_A(pm) returns an M-by-N sensing matrix according
%   to the field pm.sen_mat:
%
%     'Gaussian'         - A row of A is drawn from N(0, Sigma), with
%                          covariance Sigma = (1 - r) * I_N + r * ones(N).
%                          Required field: pm.r (default 0.5).
%
%     'Oversampled_DCT'  - A(i,j) = cos(2*pi*i*omega_j / F) / sqrt(M),
%                          with omega_j ~ U(0,1) and refinement factor F.
%                          Required field: pm.F (default 5).
%
%   This minimal implementation is sufficient to reproduce the figures
%   provided in this repository.
%
%   Reference:
%     L. Yang, S. Morigi, M. K. Ng, Y.-W. Wen,
%     "Truncated Huber Penalty for Sparse Signal Recovery with
%      Convergence Analysis," SIAM J. Sci. Comput., 48 (2026), A929-A957.

    M = pm.M;
    N = pm.N;

    switch pm.sen_mat
        case 'Gaussian'
            if isfield(pm, 'r'); r = pm.r; else; r = 0.5; end
            % Covariance: r off-diagonal, 1 on-diagonal.
            Sigma = r * ones(N) + (1 - r) * eye(N);
            A     = mvnrnd(zeros(1, N), Sigma, M);

        case 'Oversampled_DCT'
            if isfield(pm, 'F'); F = pm.F; else; F = 5; end
            w     = rand(M, 1);                          % column vector in [0,1]
            cols  = ones(M, 1) * (1:N);                  % index grid
            A     = cos(2 * pi * (w * ones(1, N)) .* cols / F) / sqrt(M);

        otherwise
            error('data_generator_A:unknownMatrix', ...
                  'Unknown pm.sen_mat: %s', pm.sen_mat);
    end
end
