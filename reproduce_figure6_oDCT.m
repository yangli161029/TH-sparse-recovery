% reproduce_figure6_DCT.m
%
% Reproduces the TH curve of Figure 6 (oversampled DCT sensing matrix) from:
%   L. Yang, S. Morigi, M. K. Ng, Y.-W. Wen,
%   "Truncated Huber Penalty for Sparse Signal Recovery with
%    Convergence Analysis," SIAM J. Sci. Comput., 48 (2026),
%    pp. A929-A957.  DOI: 10.1137/25M1748184
%
% Sensing matrix : Oversampled DCT, 64 x 1024, F = 5
% Sparsity       : s = 8
% SNR range      : 30 dB to 72 dB, step 3 dB  (15 values)
% Trials         : 50

close all; clear all;

%% ---- add path to functions folder ------------------------------------
siblingFolderName = 'functions';
currentScriptFullPath = mfilename('fullpath');
[currentScriptFolder, ~, ~] = fileparts(currentScriptFullPath);
[parentFolder, ~, ~] = fileparts(currentScriptFolder);
pathToSiblingFolder = fullfile(parentFolder, siblingFolderName);
addpath(genpath(pathToSiblingFolder));

%% ---- experimental parameters -----------------------------------------
pm.M        = 64;
pm.N        = 512 * 2;   % 1024
pm.sparsity = 8;

mindB    = 30;
maxdB    = 72;
inter    = 3;
snr_list = mindB : inter : maxdB;   % [30 33 36 ... 72], length = 15
I        = length(snr_list);
trials   = 50;

% Oversampled DCT sensing matrix settings
pm.sen_mat    = 'Oversampled_DCT';
pm.F          = 5;
pm.cond1      = 0;
pm.dynmic     = 0;
pm.normalized = 1;
pm.addon      = mindB - 3;

% TH hyper-parameters (from noisyODCTSNRsubmissionMay.m)
DataRegP_list = [2e3, 2e3, 4e3, 2e4, 6e4, 6e4, 8e4, 1e5, ...
                 5e5, 7e5, 3e6, 3e6, 8e6, 2e7, 8e6];
num_list      = [6, 3, 3, 3, 5, 8, 7, 7, 8, 7, 8, 9, 9, 9, 9];

%% ---- storage ---------------------------------------------------------
TH_error = zeros(I, trials);

%% ---- main loop -------------------------------------------------------
fprintf('Running TH on oversampled DCT matrix (64x1024), s=8, %d trials ...\n', trials);
for i = 1 : I
    dB = snr_list(i);
    fprintf('  SNR = %d dB  (%d/%d)\n', dB, i, I);

    for nt = 1 : trials
        rng(nt * 20);

        %% generate sensing matrix and sparse signal
        A     = data_generator_A(pm);
        xg    = data_generator_xg(pm);
        pm.xg = xg;
        b_gt  = A * xg;
        b     = awgn(b_gt, dB);

        %% L1 warm start
        pmL1.lambda = 1e-6;
        pmL1.maxit  = 5 * pm.N;
        pmL1.xg     = xg;
        pmL1.reltol = 1e-8;
        [x0, ~] = ADMM_L1(A, b, pmL1);

        %% TH
        pm_th.x0       = x0;
        pm_th.iter     = 5 * pm.N;
        pm_th.DataRegP = DataRegP_list(i);
        pm_th.num      = num_list(i);

        [xTH, ~] = THuberLagrange(A, b, pm_th);
        TH_error(i, nt) = norm(xTH - xg) / norm(xg);
    end
end

%% ---- compute mean RRE ------------------------------------------------
MSETH    = mean(TH_error, 2);
snr_list = snr_list(:);

%% ---- plot ------------------------------------------------------------
fontsizea = 30;
fontsizel = 40;

indices = 1 : 3 : length(snr_list);

figure;
axes('Position', [0.12, 0.18, 0.82, 0.74]);
plot(snr_list, MSETH, '-k>', ...
     'LineWidth', 3.0, 'MarkerSize', 16, 'MarkerIndices', indices);

set(gca, 'FontSize', fontsizea);
set(gca, 'YScale', 'log');
xticks(snr_list(indices));
box on;

ylabel('RRE', 'Interpreter', 'latex', 'FontSize', fontsizel, 'FontWeight', 'bold');
xlabel('SNR (dB)', 'Interpreter', 'latex', 'FontSize', fontsizel, 'FontWeight', 'bold');

legend('TH', 'Interpreter', 'latex', 'FontSize', 27, ...
       'Location', 'southwest', 'FontWeight', 'bold');

%% ---- save figure -----------------------------------------------------
print('-dpng', 'figure6_DCT_TH.png', '-r300');
fprintf('Figure saved: figure6_DCT_TH.png\n');