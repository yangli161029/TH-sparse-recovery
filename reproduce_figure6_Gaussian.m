% reproduce_figure6_Gaussian.m
%
% Reproduces the TH curve of Figure 6 (Gaussian sensing matrix) from:
%   L. Yang, S. Morigi, M. K. Ng, Y.-W. Wen,
%   "Truncated Huber Penalty for Sparse Signal Recovery with
%    Convergence Analysis," SIAM J. Sci. Comput., 48 (2026),
%    pp. A929-A957.  DOI: 10.1137/25M1748184
%
% Sensing matrix : Gaussian, 64 x 512, r = 0.5
% Sparsity       : s = 12
% SNR range      : 18 dB to 60 dB, step 3 dB  (15 values)
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
pm.N        = 512;
pm.sparsity = 12;

mindB    = 18;
maxdB    = 60;
inter    = 3;
snr_list = mindB : inter : maxdB;   % [18 21 24 ... 60], length = 15
I        = length(snr_list);
trials   = 50;

% Gaussian sensing matrix settings
pm.sen_mat    = 'Gaussian';
pm.matrixtype = 1;
pm.supptype   = 1;
pm.r          = 0.5;
pm.normalized = 0;
pm.addon      = mindB - 3;

% TH hyper-parameters (from noisyGauSNRsubmissionMay.m)
DataRegP_list = [5e1, 1e2, 2e2, 3e2, 1e3, 1e3, 3e3, 5e3, ...
                 8e3, 8e3, 2e4, 2e4, 8e3, 2e3, 2e3];
num_list      = [7, 9, 13, 10, 12, 12, 9, 8, 13, 9, 9, 12, 13, 13, 13];

%% ---- storage ---------------------------------------------------------
TH_error = zeros(I, trials);

%% ---- main loop -------------------------------------------------------
fprintf('Running TH on Gaussian matrix (64x512), s=12, %d trials ...\n', trials);
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

indices = 1 : 2 : length(snr_list);

figure;
axes('Position', [0.12, 0.18, 0.82, 0.74]);
plot(snr_list, MSETH, '-k>', ...
     'LineWidth', 3.0, 'MarkerSize', 16, 'MarkerIndices', indices);

set(gca, 'FontSize', fontsizea);
set(gca, 'YScale', 'log');
xticks(snr_list);
box on;

ylabel('RRE', 'Interpreter', 'latex', 'FontSize', fontsizel, 'FontWeight', 'bold');
xlabel('SNR (dB)', 'Interpreter', 'latex', 'FontSize', fontsizel, 'FontWeight', 'bold');

legend('TH', 'Interpreter', 'latex', 'FontSize', 27, ...
       'Location', 'southwest', 'FontWeight', 'bold');

%% ---- save figure -----------------------------------------------------

print('-dpng',  'figure6_Gaussian_TH.png', '-r300');
fprintf('Figure saved: figure6_Gaussian_TH.{png}\n');
