%REPRODUCE_FIGURE3_DCT  Reproduce the TH panel of Fig. 3 (Oversampled DCT).
%
%  Paper:  L. Yang, S. Morigi, M. K. Ng, Y.-W. Wen,
%          "Truncated Huber Penalty for Sparse Signal Recovery with
%           Convergence Analysis," SIAM J. Sci. Comput., 48 (2026),
%           pp. A929-A957.    DOI: 10.1137/25M1748184
%
%  This script reproduces the ground-truth panel and the TH panel of
%  the bottom three rows of Figure 3.  The original figure additionally
%  shows seven baseline methods; those baseline implementations are
%  NOT redistributed here -- please obtain them from their respective
%  authors / public packages.
%
%  Setting (noise-free):
%     A2 in R^{64 x 1024},  oversampled DCT with F = 5
%     ground truth: 14 nonzero entries with exponential decay,
%                   placed at random positions.
%
%  Expected output:
%     RRE          ~ 6.6e-16
%     iterations   ~ 34
%     runtime      < 0.1 s on a recent laptop
%     PNG saved to ./outputs/figure3_TH_DCT.png

close all; clear; clc;

%% --- Path setup (relative) -----------------------------------------
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
addpath(genpath(fullfile(thisDir, 'functions')));

outputDir = fullfile(thisDir, 'outputs');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

%% --- Problem setup --------------------------------------------------
pm.M        = 64;
pm.N        = 1024;
pm.tau      = 0;                 % noise-free
pm.sparsity = 14;
pm.sen_mat  = 'Oversampled_DCT';
pm.F        = 5;
pm.restol   = 1e-8;

rng(4*20);                       % reproducibility (matches paper run)
A = data_generator_A(pm);

% Ground truth: exponentially decaying, then permuted
xg          = zeros(pm.N, 1);
idx         = 1:15:200;
xg(idx)     = sqrt(200) ./ (1 + exp((idx - 100) / 20));
xg          = xg(randperm(pm.N));
pm.xg       = xg;

% Noise-free measurement
b = A * xg;

%% --- L1 ADMM warm start --------------------------------------------
pmL1.lambda = 1e-6;
pmL1.maxit  = 5 * pm.N;
pmL1.xg     = xg;
pmL1.reltol = 1e-8;
[xL1, ~]    = ADMM_L1(A, b, pmL1);

%% --- TH (truncated Huber) recovery ---------------------------------
pm.x0       = xL1;
pm.iter     = 5 * pm.N;
pm.DataRegP = 3e3;               % continuation: fidelity weight (DCT)
pm.num      = 13;                % continuation: number of epochs / index
pm.restol   = 1e-8;

tStart       = tic;
[xTH, outTH] = THuberLagrange(A, b, pm);
xTHtime      = toc(tStart);
xTHitecount  = outTH.iter_count;
xTHrre       = norm(xTH - xg) / norm(xg);

fprintf('=== Figure 3 (Oversampled DCT) -- TH ===\n');
fprintf('  RRE        : %.4e\n', xTHrre);
fprintf('  Runtime (s): %.4f\n', xTHtime);
fprintf('  Iterations : %d\n',   xTHitecount);

%% --- Plot (1x2 panel: ground truth + TH) ---------------------------
fontsizet = 22;
fontsizea = 16;
titlePos  = [0.5, -0.18, 0];

fig = figure('Position', [100, 100, 1200, 420]);

% (a) ground truth
ax = axes('Position', [0.06, 0.22, 0.42, 0.60]); %#ok<LAXES>
stem(max(abs(xg), 1e-10), 'MarkerSize', 8, 'LineWidth', 1.4);
set(ax, 'YScale', 'log', 'FontSize', fontsizea, ...
        'XLim', [1, pm.N], 'YLim', [1e-10, 1e2]);
text(0.04, 0.92, '(a)', 'Units', 'normalized', ...
     'FontSize', fontsizea, 'FontWeight', 'bold');
t = title(ax, '\textbf{ground truth}', 'Interpreter', 'latex', ...
          'FontSize', fontsizet);
set(t, 'Units', 'normalized', 'Position', titlePos, ...
       'VerticalAlignment', 'top');

% (b) TH reconstruction
ax = axes('Position', [0.55, 0.22, 0.42, 0.60]); %#ok<LAXES>
stem(max(abs(xTH), 1e-10), 'MarkerSize', 8, 'LineWidth', 1.4);
set(ax, 'YScale', 'log', 'FontSize', fontsizea, ...
        'XLim', [1, pm.N], 'YLim', [1e-10, 1e2]);
text(0.04, 0.92, '(b)', 'Units', 'normalized', ...
     'FontSize', fontsizea, 'FontWeight', 'bold');
t = title(ax, sprintf('\\textbf{TH}: %.2e, %.2fs, %d', ...
                      xTHrre, xTHtime, xTHitecount), ...
          'Interpreter', 'latex', 'FontSize', fontsizet);
set(t, 'Units', 'normalized', 'Position', titlePos, ...
       'VerticalAlignment', 'top');

% Save
outFile = fullfile(outputDir, 'figure3_TH_DCT.png');
if exist('exportgraphics', 'file') == 2
    exportgraphics(fig, outFile, 'Resolution', 300);
else
    print(fig, outFile, '-dpng', '-r300');
end
fprintf('\nFigure saved to: %s\n', outFile);
