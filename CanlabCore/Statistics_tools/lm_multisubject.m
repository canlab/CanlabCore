function [coeff_matrix, Rsq, d_effect_size_matrix, Rsq_partial] = lm_multisubject(model, varargin)
% lm_multisubject performs subject-level linear regression analyses and
% creates visualizations for fixed-effects coefficients, R-squared values, 
% partial R-squared values, and effect sizes across multiple subjects.
%
% Inputs:
%   - model: A cell array of fitted models (from `fitlm`), one for each subject.
%   - varargin (optional, but both must be provided):
%       - sub_regressors: A cell array containing the predictor variables used for 
%         each subject's model.
%       - sub_y: A cell array containing the outcome variables (Y values) for each subject.
%
% Outputs:
%   - coeff_matrix: A matrix of coefficients for each subject and regressor.
%   - Rsq: A matrix containing the Ordinary and Adjusted R-squared values for each subject.
%   If varargin is provided:
%   - d_effect_size_matrix: A matrix of Cohen's d effect sizes for each subject and regressor.
%   - Rsq_partial: A matrix of partial R-squared values for each subject and regressor.
%
% Visualizations:
%   - Heatmaps of subject-level fixed-effects coefficients with significance indicators.
%   - Floating error bars for each subject's coefficients.
%   - Heatmaps of subject-level R-squared values (Ordinary and Adjusted).
%   If varargin is provided:
%   - Heatmaps of subject-level effect sizes (Cohen's d) and partial R-squared values.
%
% Example usage:
%   [coeff_matrix, Rsq, d_effect_size_matrix, Rsq_partial] = lm_multisubject(models, sub_regressors, sub_y);
%
% The function automatically generates plots that visualize the following:
%   - Fixed-effects coefficients for each subject (with significance indicators).
%   - Floating error bars for each subject's coefficients.
%   - R-squared heatmaps showing the model fit across subjects.
%   - Partial R-squared values and effect sizes (Cohen's d) for each regressor and subject.
%
% Notes:
%   - The function can take in both `sub_regressors` and `sub_y` to compute partial R-squared
%     and effect size matrices.
%
% Author: Michael Sun, Ph.D. 08/15/2024

d_effect_size_matrix=[];
Rsq_partial=[];
noplot=0;

% Handle optional inputs
if ~isempty(varargin)
    optionalargs=1;
    sub_regressors = varargin{1};
    sub_y = varargin{2};
else
    optionalargs=0;
end

if any(strcmp(varargin, 'noplot'))
    noplot=1;
end

% Initialize Variables
num_subjects = numel(model);
regnames = model{1}.CoefficientNames(2:end); % Exclude intercept
num_regressors = numel(regnames);

% Pre-allocate matrices
coeff_matrix = zeros(num_subjects, num_regressors);
p_value_matrix = zeros(num_subjects, num_regressors);
conf_int_matrix = zeros(num_subjects, num_regressors, 2); % Confidence intervals
Rsq = zeros(2, num_subjects); % R2 (Ordinary and Adjusted)
Rsq_partial = zeros(num_subjects, num_regressors);
d_effect_size_matrix = zeros(num_subjects, num_regressors);

% Extract coefficients, p-values, and R-squared values from each model
for i = 1:num_subjects
    coeff_matrix(i, :) = model{i}.Coefficients.Estimate(2:end)';
    p_value_matrix(i, :) = model{i}.Coefficients.pValue(2:end)';
    ci = coefCI(model{i});
    conf_int_matrix(i, :, 1) = ci(2:end, 1)';
    conf_int_matrix(i, :, 2) = ci(2:end, 2)';
    Rsq(1, i) = model{i}.Rsquared.Ordinary;
    Rsq(2, i) = model{i}.Rsquared.Adjusted;
end


if noplot==0
% Create Plots and Tables
figure;

% Coefficient Heatmap

if optionalargs==0
    subplot(2,2,1);
else
    subplot(3,2,1);
end

imagesc(coeff_matrix);
bound=max(abs([max(coeff_matrix), min(coeff_matrix)]));
clim([-1*bound bound])

cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]); 
colormap(cm);
colorbar;
title('Subject-Level Fixed-Effects Coefficients');
xlabel('Regressors');
ylabel('Subjects');
set(gca, 'XTick', 1:num_regressors, 'XTickLabel', regnames);
set(gca, 'YTick', 1:num_subjects, 'YTickLabel', 1:num_subjects);

% Overlay significance indicators
hold on;
for i = 1:num_subjects
    for j = 1:num_regressors
        if p_value_matrix(i, j) < 0.05
            rectangle('Position', [j-0.5, i-0.5, 1, 1], 'EdgeColor', 'black', 'LineWidth', 2);
        end
    end
end
hold off;

% Floating Error Bars
subplot(3,2,2);
hold on;
jitterAmount = 0.1;
for i = 1:num_subjects
    jitter = (rand(1, num_regressors) - 0.5) * 2 * jitterAmount;
    x_positions = (1:num_regressors) + jitter;
    errorbar(x_positions, coeff_matrix(i, :), ...
             coeff_matrix(i, :) - conf_int_matrix(i, :, 1), ...
             conf_int_matrix(i, :, 2) - coeff_matrix(i, :), ...
             'o', 'LineWidth', 1.5);
end
title('Floating Error Bars of Subject-Level Fixed Effects');
xlabel('Regressors');
ylabel('Coefficient Estimate');
xticks(1:num_regressors);
xticklabels(regnames);
hold off;

% R-squared Heatmap
subplot(3,2,3:4);
imagesc(Rsq);
bound=max(abs([max(coeff_matrix), min(coeff_matrix)]));
clim([-1*bound bound])
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]); 
colormap(cm);
colorbar;
title('Subject-Level Model R2');
xlabel('Subjects');
set(gca, 'YTick', 1:2, 'YTickLabel', {'R2 (Ord)', 'R2 (Adj)'});

% Partial R2 and Effect Size Computation
for i = 1:num_subjects
    % Full R2
    Rsq_full = model{i}.Rsquared.Ordinary;
    for j = 1:num_regressors
        % Create reduced model by excluding the j-th regressor
        reduced_regressors = [1:j-1, j+1:num_regressors];
        mdl_reduced = fitlm(sub_regressors{i}(:, reduced_regressors), sub_y{i});
        Rsq_partial(i, j) = Rsq_full - mdl_reduced.Rsquared.Ordinary;
        betas = model{i}.Coefficients.Estimate(2:end)';
        std_predictors = std(sub_regressors{i}(:, j));
        std_response = std(sub_y{i});
        d_effect_size_matrix(i, j) = (betas(j) * std_predictors) / std_response;
    end
end

if optionalargs==1

% Display Partial R2 Heatmap and Effect Sizes
subplot(3,2,5);
imagesc(d_effect_size_matrix);
bound=max(abs([max(coeff_matrix), min(coeff_matrix)]));
clim([-1*bound bound])
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]); 
colormap(cm);
colorbar;
% Overlay significance indicators
hold on;
for i = 1:num_subjects
    for j = 1:num_regressors
        if p_value_matrix(i, j) < 0.05
            rectangle('Position', [j-0.5, i-0.5, 1, 1], 'EdgeColor', 'black', 'LineWidth', 2);
        end
    end
end
title('Subject-Level Effect Sizes');
xlabel('Subjects');
ylabel('Regressors');

subplot(3,2,6);
imagesc(Rsq_partial);
bound=max(abs([max(coeff_matrix), min(coeff_matrix)]));
clim([-1*bound bound])
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]); 
colormap(cm);
colorbar;
% Overlay significance indicators
hold on;
for i = 1:num_subjects
    for j = 1:num_regressors
        if p_value_matrix(i, j) < 0.05
            rectangle('Position', [j-0.5, i-0.5, 1, 1], 'EdgeColor', 'black', 'LineWidth', 2);
        end
    end
end
title('Partial R2');
xlabel('Subjects');
ylabel('Regressors');
end
end

end