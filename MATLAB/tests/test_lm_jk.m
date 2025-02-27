%% Step 1: Simulate Data
clc; clear; close all;

% Set relative paths based on current script location
scriptDir = fileparts(mfilename('fullpath'));  % Get current script's directory
projectRoot = fullfile(scriptDir, '../..');   % Go up two levels to reach `/FMR/`
addpath(fullfile(projectRoot, 'otherfunctions')); % Add `otherfunctions` folder

% Set output directory (test folder)
outputDir = fullfile(scriptDir, 'results');  % Save results in `test/`
if ~exist(outputDir, 'dir')
    mkdir(outputDir);  % Create `test/` folder if it doesn't exist
end

% Set seed for reproducibility
rng(42);

% Define dimensions
nn = 100;  % Number of observations (SNPs)
pp = 5;    % Number of predictors (features)

% Generate predictor matrix X (standard normal)
X = randn(nn, pp);

% Define true beta coefficients
true_beta = [1.5; -2; 3; 0; 0.5];

% Generate response variable y = X * true_beta + noise
noise = 0.5 * randn(nn, 1);  % Gaussian noise (sigma = 0.5)
y = X * true_beta + noise;

% Generate random weights W (between 0.5 and 1.5)
W = 0.5 + rand(nn, 1);

%% Step 2: Apply Jackknife Regression
disp('Running Jackknife Regression...');

% Call the function
[beta_jk, exitflag] = linear_regression_jk(X, y, 'W', W, 'NoJackknifeBlocks', 10);

% Compute mean Jackknife estimate
beta_mean = mean(beta_jk, 1);

%% Step 3: Compare with Standard Regression
disp('Running Standard Weighted Least Squares Regression...');

% Standard weighted regression solution
beta_standard = (X' * diag(W) * X) \ (X' * diag(W) * y);

%% Step 4: Display Results
disp('True Beta Coefficients:');
disp(true_beta');

disp('Standard Weighted Least Squares Estimate:');
disp(beta_standard');

disp('Mean Jackknife Regression Estimate:');
disp(beta_mean);

% Compute standard error of Jackknife estimates
jk_std_err = std(beta_jk, 0, 1);

disp('Jackknife Standard Error:');
disp(jk_std_err);

%% Step 5: Save Results
save(fullfile(outputDir, 'jackknife_results.mat'), 'beta_mean', 'beta_jk', 'true_beta', 'jk_std_err');

disp(['Results saved in: ', fullfile(outputDir, 'jackknife_results.mat')]);

%% Step 6: Plot and Save Figures
figure;
hold on;
errorbar(1:pp, beta_mean, jk_std_err, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(1:pp, true_beta, 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
xlabel('Coefficient Index');
ylabel('Estimate');
legend('Jackknife Mean ± SE', 'True Beta', 'Location', 'best');
title('Jackknife Regression Estimates');
grid on;
hold off;

% Save figure
saveas(gcf, fullfile(outputDir, 'jackknife_plot.png'));

disp(['Figure saved in: ', fullfile(outputDir, 'jackknife_plot.png')]);
