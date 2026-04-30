% reproduce_figure8.m
%
% Reproduces the TH-related panels of Figure 8 from:
%   L. Yang, S. Morigi, M. K. Ng, Y.-W. Wen,
%   "Truncated Huber Penalty for Sparse Signal Recovery with
%    Convergence Analysis," SIAM J. Sci. Comput., 48 (2026),
%    pp. A929-A957.  DOI: 10.1137/25M1748184
%
% Panels reproduced:
%   (a) Noisy input signal
%   (d) Recovered auxiliary variable omega (w)
%   (f) TH reconstruction result
%
% Signal  : 1D piecewise constant (blocks), N = 256
% Noise   : additive white Gaussian noise (sigma = 0.5)
% Data    : data/blocks.txt  (clean signal)
%           data/blocks_noisy.txt  (noisy signal)
%
% Required functions (all in the functions/ folder):
%   algvetgradiffbdy.m, Diffmatrixboundarycond.m, vecsubprb.m

close all; clear all;

%% ---- add path to functions folder ------------------------------------
siblingFolderName = 'functions';
currentScriptFullPath = mfilename('fullpath');
[currentScriptFolder, ~, ~] = fileparts(currentScriptFullPath);
pathToSiblingFolder = fullfile(currentScriptFolder, siblingFolderName);
addpath(genpath(pathToSiblingFolder));

%% ---- create output folder --------------------------------------------
outputFolder = fullfile(currentScriptFolder, 'output');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% ---- load data -------------------------------------------------------
dataFolder = fullfile(currentScriptFolder, 'data');
x_clean    = load(fullfile(dataFolder, 'blocks.txt'));       % clean signal
y          = load(fullfile(dataFolder, 'blocks_noisy.txt')); % noisy signal
N          = length(y);
n          = 1 : N;

%% ---- TH parameters (from denoiseforsubmission.m) ---------------------
lambda = 1e4;
thr    = 1.5e0;
maxit  = 200;
flag   = 4;        % Dirichlet boundary condition

%% ---- run TH ----------------------------------------------------------
fprintf('Running TH denoising (N=%d) ...\n', N);
[unew, w, result] = algvetgradiffbdy(y, lambda, thr, maxit, flag);
fprintf('Done. Outer iterations: %d\n', length(result.i_list));

%% ---- metrics ---------------------------------------------------------
err_rmse = sqrt(mean((x_clean - unew).^2));
err_mae  = mean(abs(x_clean - unew));
fprintf('RMSE: %.4f\n', err_rmse);
fprintf('MAE : %.4f\n', err_mae);

%% ---- plot settings ---------------------------------------------------
ax        = [0 256 -3 6];
fontsize  = 25;
linewidth = 3.0;

%% ---- panel (a): noisy signal -----------------------------------------
figure;
axes('Position', [0.12, 0.18, 0.82, 0.74]);
plot(n, y, 'k-', 'LineWidth', linewidth);
axis(ax);
set(gca, 'FontSize', fontsize, 'LineWidth', linewidth);
box on;
xlabel('Index', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', fontsize);
title('(a) Noisy signal', 'Interpreter', 'latex', 'FontSize', fontsize);
savePath = fullfile(outputFolder, 'figure8a_noisy.png');
print('-dpng', savePath, '-r300');
fprintf('Figure saved: %s\n', savePath);

%% ---- panel (d): recovered omega (w) ----------------------------------
figure;
axes('Position', [0.12, 0.18, 0.82, 0.74]);
stem(1 - w, ...
    'MarkerFaceColor', 'b', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerSize', 6, ...
    'LineStyle', 'none');
axis([0 256 -0.2 1.2]);
set(gca, 'FontSize', fontsize);
grid on;
box on;
xlabel('Index', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('$\hat{\omega}$', 'Interpreter', 'latex', 'FontSize', fontsize);
title('(d) Recovered $\hat{\omega}$', 'Interpreter', 'latex', 'FontSize', fontsize);
savePath = fullfile(outputFolder, 'figure8d_omega.png');
print('-dpng', savePath, '-r300');
fprintf('Figure saved: %s\n', savePath);

%% ---- panel (f): TH reconstruction ------------------------------------
figure;
axes('Position', [0.12, 0.18, 0.82, 0.74]);
plot(n, x_clean, '--', 'LineWidth', linewidth);
hold on;
plot(n, unew,    '-',  'LineWidth', linewidth);
hold off;
axis(ax);
set(gca, 'FontSize', fontsize, 'LineWidth', linewidth);
box on;
xlabel('Index', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', fontsize);
title('(f) TH', 'Interpreter', 'latex', 'FontSize', fontsize);
text(0.05, 0.9, sprintf('RMSE: %.4f\nMAE: %.4f', err_rmse, err_mae), ...
    'Units', 'normalized', 'FontSize', fontsize, 'Interpreter', 'latex');
legend({'Ground truth', 'TH'}, 'Interpreter', 'latex', ...
    'FontSize', 20, 'Location', 'northeast');
savePath = fullfile(outputFolder, 'figure8f_TH.png');
print('-dpng', savePath, '-r300');
fprintf('Figure saved: %s\n', savePath);