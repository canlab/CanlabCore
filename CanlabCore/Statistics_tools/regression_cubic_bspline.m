function out = regression_cubic_bspline(X, y, varargin)
% REGRESSION_CUBIC_BSPLINE - Fit a cubic B-spline regression model with column expansion.
%
% :Usage:
% ::
%
%    out = regression_cubic_bspline(X, y, [optional inputs])
%
% :Inputs:
%
%    **X** : [NxM double]
%         Numeric matrix of predictor variables.
%           or
%         table object with numeric predictor variables.
%         column names will be used as predictor names.
%
%    **y** : [Nx1 double]
%         Vector of outcome values.
%
% :Optional Inputs:
%
%    **'expand_cols'** : [logical scalar or 1xM logical vector]
%         Specifies which columns of X to expand with a cubic B-spline basis.
%         If a scalar and true, all columns are expanded; if false, none are expanded.
%         *Default:* true.
%
%    **'names'** : [cell array of strings]
%         Names for each column of X.
%         *Default:* {} (if empty, columns are named 'Column_1', 'Column_2', etc.).
%
%    **'verbose'** : [logical scalar]
%         Flag to print progress information.
%         *Default:* true.
%
%    **'doplot'** : [logical scalar]
%         Flag to produce plots (actual vs. predicted, partial regression plots,
%         and basis function plots).
%         *Default:* true.
%
% :Outputs:
%
%    **out** : [structure]
%         Structure with the following fields:
%
%         - **X**: [NxP double] Final design matrix (including intercept as the last column).
%         - **beta**: [Px1 double] Fitted regression coefficients.
%         - **y_hat**: [Nx1 double] Predicted responses.
%         - **partial_residuals**: [1xM cell] Partial residuals for each original regressor.
%         - **F_tests**: [table] F-statistics for each original regressor, with columns:
%               Regessor, F_stat, df1, df2, and p_value.
%         - **knot_pts**: [structure] Knot points used for each expanded column of X.
%         - **basis_sets**: [structure] B-spline basis matrices for each expanded column of X.
%         - **expanded_regs**: [vector] Indices of regressors (columns of X) that were expanded.
%         - **mapping_matrix**: [logical matrix]
%               Mapping from original regressors to columns of the expanded design matrix
%               (before adding the intercept).
%         - **input_params**: [structure] All input parameters passed to the function.
%         - **names**: [cell array of strings] Names of regressors in the final design matrix.
%
% :Examples:
%
%    % Generate data for age, practice, and memory:
%    rng(0);          % For reproducibility
%    N = 100;         % Sample size
%
%    % Generate 'age' uniformly between 20 and 80:
%    age = 20 + 60 * rand(N, 1);
%
%    % Generate 'practice' uniformly between 5 and 200:
%    practice = 5 + (200 - 5) * rand(N, 1);
%
%    % Generate 'memory' as a nonlinear function of age and practice.
%    % Memory is modeled as a quadratic function of age and practice with added noise.
%    memory = -0.02 * (age - 30).^2 + 0.002 * (practice - 100).^2 + 50 + randn(N, 1) * 3;
%
%    % Check that age and practice are not correlated:
%    corr_age_practice = corr(age, practice);
%    fprintf('Linear correlation between age and practice: %.3f\n', corr_age_practice);
%
%    % Example 1: Single regressor ("age") expanded:
%    out = regression_cubic_bspline(age, memory, 'expand_cols', true, 'names', {'age'}, ...
%          'verbose', true, 'doplot', true);
%
%    % Example 2: Two regressors ("age" and "practice hours")
%
%    % Case A: Only age is expanded:
%    out1 = regression_cubic_bspline([age practice], memory, 'expand_cols', [true, false], ...
%           'names', {'age', 'practice_hours'}, 'verbose', true, 'doplot', true);
%
%    % Case B: Both age and practice hours are expanded:
%    out2 = regression_cubic_bspline([age practice], memory, 'expand_cols', [true, true], ...
%           'names', {'age', 'practice_hours'}, 'verbose', true, 'doplot', true);
%
% :References:
%
%    Uses some functions from Canlab Core Tools
%    (https://github.com/canlab)
%    And Matlab curve fitting toolbox
%
% :See also:
%
%    augknt, spcol in curve fitting toolbox for key spline functions
%
% Author: Tor Wager

% Copyright (C) 2025  Tor Wager
% GNU General Public License; see http://www.gnu.org/licenses/

%% Input Parsing and Validation
p = inputParser;

addRequired(p, 'X', @(z) validateattributes(z, {'table' 'numeric'}, {'2d', 'nonempty'}));

addRequired(p, 'y', @(z) validateattributes(z, {'numeric'}, {'column'}));

addParameter(p, 'expand_cols', true, @(z) islogical(z) || isnumeric(z));

addParameter(p, 'names', {}, @(z) iscell(z) && all(cellfun(@ischar, z)));

addParameter(p, 'verbose', true, @(z) islogical(z) && isscalar(z));

addParameter(p, 'doplot', true, @(z) islogical(z) && isscalar(z));

parse(p, X, y, varargin{:});

input_params = p.Results;

[n_rows, n_cols] = size(X);

if istable(X)
    input_params.names = X.Properties.VariableNames;
    X = table2array(X);
end

validateattributes(y, {'numeric'}, {'numel', n_rows});


% Process expand_cols input: if scalar, expand all (true) or none (false)
if islogical(input_params.expand_cols) && isscalar(input_params.expand_cols)

    expand_cols = repmat(input_params.expand_cols, 1, n_cols);

elseif islogical(input_params.expand_cols) && numel(input_params.expand_cols) == n_cols

    expand_cols = input_params.expand_cols;

else
    error('expand_cols must be a scalar logical or a logical vector with length equal to the number of columns in X.');
end


% Process names input: if empty, create default names.
if ~isempty(input_params.names)

    if numel(input_params.names) ~= n_cols
        error('The length of names must equal the number of columns in X.');
    end

    names = input_params.names;

else
    names = cell(1, n_cols);
    for i = 1 : n_cols
        names{i} = sprintf('Column_%d', i);
    end
end


%% Expand Selected Columns using Cubic B-Spline Basis
knot_pts      = struct();
basis_sets    = struct();
expanded_regs = [];

new_columns   = {};
new_names     = {};
orig_idx      = cell(1, n_cols);

col_counter   = 0;

for i = 1 : n_cols

    col_data = X(:, i);

    if input_params.doplot

        fig_han(i) = create_figure(sprintf('Plots for %s', names{i}), 1, 2);

    end

    if expand_cols(i)

        if input_params.verbose
            fprintf('Expanding column %d (%s) with cubic B-spline basis...\n', i, names{i});
        end

        [basis, knots] = expand_column_with_spline(col_data, input_params.verbose);

        knot_pts.(sprintf('col_%d', i))   = knots;
        basis_sets.(sprintf('col_%d', i)) = basis;

        expanded_regs = [expanded_regs, i];

        n_basis = size(basis, 2);
        indices = col_counter + (1 : n_basis);
        col_counter = col_counter + n_basis;
        orig_idx{i} = indices;

        for j = 1 : n_basis
            new_columns{end+1} = basis(:, j);
            new_names{end+1}   = sprintf('%s_bspline_%d', names{i}, j);
        end

        % Plot the basis functions for the expanded column
        if input_params.doplot

            [col_data_sorted, sort_idx] = sort(col_data);
            basis_sorted = basis(sort_idx, :);

            for j = 1 : n_basis
                plot(col_data_sorted, basis_sorted(:, j), 'LineWidth', 1.5)
            end

            hold off
            xlabel(names{i});
            ylabel('Basis Value');
            title(sprintf('Basis Functions for %s', names{i}));
        end

    else

        if input_params.doplot

            plot(col_data, col_data, 'LineWidth', 1.5)

        end

        indices = col_counter + 1;
        col_counter = col_counter + 1;
        orig_idx{i} = indices;

        new_columns{end+1} = col_data;
        new_names{end+1}   = names{i};

    end

end


% Create a mapping matrix (logical) for the non-intercept regressors.
mapping_matrix = false(n_cols, col_counter);
for i = 1 : n_cols
    mapping_matrix(i, orig_idx{i}) = true;
end


% Assemble the final design matrix using new_columns
X_expanded = [new_columns{:}];


%% Always add an intercept as the last column
X_expanded = intercept(X_expanded, 'end');

new_names{end + 1} = 'Intercept';

%% Fit Linear Regression Model
beta  = X_expanded \ y;
y_hat = X_expanded * beta;


%% Partial Regression Plots and F-tests for Each Original Regressor
residuals = y - y_hat;

partial_residuals = cell(1, n_cols);

F_reg_names = cell(n_cols, 1);
F_stats     = zeros(n_cols, 1);
df1s        = zeros(n_cols, 1);
df2s        = zeros(n_cols, 1);
p_values    = zeros(n_cols, 1);

for i = 1 : n_cols


    % Compute F-test by comparing full model to reduced model (excluding regressor i)
    % ---------------------------------------------------------
    X_reduced = X_expanded;
    X_reduced(:, orig_idx{i}) = [];

    beta_reduced = X_reduced \ y;
    RSS_reduced = sum((y - X_reduced * beta_reduced).^2);
    RSS_full    = sum(residuals.^2);

    p_i = numel(orig_idx{i});
    n_obs = n_rows;
    p_full = size(X_expanded, 2);

    df1 = p_i;
    df2 = n_obs - p_full;

    F_stat = ((RSS_reduced - RSS_full) / df1) / (RSS_full / df2);

    F_reg_names{i} = names{i};
    F_stats(i) = F_stat;
    df1s(i) = df1;
    df2s(i) = df2;
    p_values(i) = 1 - fcdf(F_stat, df1, df2);

    if input_params.verbose
        fprintf('F-test for regressor %d (%s): F = %.3f, df1 = %d, df2 = %d, p = %.4f\n', ...
            i, names{i}, F_stat, df1, df2, p_values(i));
    end

    % Compute variables for partial residual plots and output
    % ---------------------------------------------------------
    % Compute the fitted effect for the original regressor using its associated columns
    effect_i = X_expanded(:, orig_idx{i}) * beta(orig_idx{i});

    % The fit line for the partial regression plot is the sum of the regressor's effect
    % and the intercept (last column of beta)
    fit_line = effect_i + beta(end);

    beta_adj_pred = X_reduced \ X(:, i);
    adj_pred{i} = X(:, i) - X_reduced * beta_adj_pred + beta_adj_pred(end); % adjusted predictor for other covariates, adding in the intercept

    partial_residuals{i} = effect_i + residuals + beta(end);


    if input_params.doplot
        % For partial regression, plot the original regressor vs. the outcome.
        % Compute the fit line using all associated expanded columns times beta, plus the intercept.
        orig_regressor = X(:, i);
        adj_regressor = adj_pred{i};

        % Sort for plotting the fit line
        [x_sorted, sort_idx] = sort(orig_regressor);
        % adj_regressor_sorted = adj_regressor(sort_idx);
        fit_line_sorted = fit_line(sort_idx);

        figure(fig_han(i))
        subplot(1, 2, 2)
        % h_partial = create_figure(sprintf('Partial Regression Plot for %s', names{i}));
        hold on

        % scatter(orig_regressor, y, 'filled')
        scatter(adj_regressor, partial_residuals{i}, 'filled')
        plot(x_sorted, fit_line_sorted, 'r-', 'LineWidth', 1.5)

        hold off
        xlabel([names{i} ' (Adjusted)']);
        ylabel('Outcome (Adjusted)');
        title(sprintf('Partial Regression Plot for %s', names{i}));
    end

end

F_tests = table(F_reg_names, F_stats, df1s, df2s, p_values, ...
    'VariableNames', {'Regressor', 'F_stat', 'df1', 'df2', 'p_value'});


%% Plot Actual vs. Predicted if Requested
if input_params.doplot

    h_actual = create_figure('Design matrix and fits', 1, 2);

    subplot(1, 2, 1)
    X_to_plot = X_expanded ./ max(abs(X_expanded)); % normalize for plot
    imagesc(X_to_plot); set(gca, 'YDir', 'reverse'); axis tight
    colormap(colormap_tor([0 0 1], [1 1 0], [.2 .5 1], [.5 .5 .5], [1 .5 .2]))
    title('Design matrix with expanded regressors')

    names_to_plot = strrep(new_names, '_b', ' b');
    set(gca, 'XTick', 1:length(new_names), 'XTickLabelRotation', 45, 'XtickLabel', names_to_plot)
    ylabel('Observation')

    subplot(1, 2, 2)
    hold on

    scatter(y_hat, y, 'filled')
    lims = [min(min(y_hat), min(y)), max(max(y_hat), max(y))];
    plot(lims, lims, 'r--', 'LineWidth', 1.5)

    hold off
    xlabel('Predicted y\_hat')
    ylabel('Actual y')
    title('Actual vs. Predicted (Cubic B-spline Regression)')

end


%% Assemble Output Structure
out = struct();

out.X                 = X_expanded;
out.beta              = beta;
out.y_hat             = y_hat;
out.partial_residuals = partial_residuals;
out.F_tests           = F_tests;
out.knot_pts          = knot_pts;
out.basis_sets        = basis_sets;
out.expanded_regs     = expanded_regs;
out.mapping_matrix    = mapping_matrix;
out.input_params      = input_params;
out.names             = new_names;


if input_params.verbose

    fprintf('Cubic B-spline regression completed. Model fitted with %d regressors.\n', size(X_expanded, 2));

    disp(out.F_tests)
end

end


% =========================================================================
function [basis, knots] = expand_column_with_spline(x_col, verbose)
% EXPAND_COLUMN_WITH_SPLINE Expand a predictor column using a cubic B-spline basis.
%
% USAGE:
%    [basis, knots] = expand_column_with_spline(x_col, verbose)
%
% INPUTS:
%    x_col   - [Nx1 double] Vector of predictor values.
%    verbose - [logical scalar] If true, prints information about knot selection.
%
% OUTPUTS:
%    basis   - [NxK double] Matrix of cubic B-spline basis functions evaluated at x_col.
%    knots   - [vector] The augmented knot vector used for constructing the spline basis.
%
% DESCRIPTION:
%    The function automatically chooses two interior knots at the 33rd and 66th
%    percentiles of x_col, constructs an augmented knot vector for a cubic spline
%    (order 4), and then computes the B-spline basis using spcol.
%
% EXAMPLE:
%    x = linspace(0, 10, 100)';
%    [basis, knots] = expand_column_with_spline(x, true);
%
% SEE ALSO:
%    augknt, spcol
%

x_col = x_col(:);

N = length(x_col);

% When your data are skewed, it is often preferable to choose knot points 
% based on percentiles of the data rather than using linear spacing across 
% the range. Using percentiles tailors the knot placement to the actual 
% density of the data, ensuring that there are more knots where the data 
% are concentrated and fewer where there are sparse observations. This 
% approach can lead to a better fit in regions with high data density and 
% helps avoid overfitting in areas with few observations. 

k1 = prctile(x_col, 33);
k2 = prctile(x_col, 66);

if verbose
    fprintf('    Chosen interior knots at %.3f and %.3f.\n', k1, k2);
end

% In standard practice, the boundary knots of a B‑spline basis are set to 
% the minimum and maximum values of the data to ensure the spline “covers” 
% the range of the data. Typically, these boundary knots are then augmented 
% (i.e., repeated a number of times equal to the spline order) to produce a 
% clamped spline that has desirable boundary properties (such as 
% interpolating the endpoints).
% 
% However, there are cases where one might choose knots that extend slightly beyond the data range, such as when the analyst wishes to avoid boundary effects or when known characteristics of the underlying process suggest that the data are only a subset of the full range of variability. Extending the knot sequence beyond the range of the data can improve extrapolation or stabilize the spline near the edges, but it also may introduce artifacts if the extrapolated behavior is not justified by the data.

knots_vec = [min(x_col), k1, k2, max(x_col)];
order = 4;
knots = augknt(knots_vec, order);

% ties in data will break the construction of the spline collocation vectors (regressors)
% if there are no unique knot points. 
% if the data are not truly continuous and there are mulitple exact repeats of data values
% this will break the spline basis construction. add a bit of noise in this
% case
if length(unique(x_col)) < size(x_col, 1)
    x_col = x_col + 0.001 * range(x_col) * unifrnd(-1, 1, size(x_col, 1), 1);
end

% we could also use unique values instead, and reconstruct X_sorted

[x_sorted, sort_idx] = sort(x_col);
X_sorted = spcol(knots, order, x_sorted);  % Create B-spline collocation basis

basis = zeros(N, size(X_sorted, 2));
basis(sort_idx, :) = X_sorted;

end

% Notes:
% repeated values (ties) will break the spline construction
% Error using spcol (line 121)
% Point multiplicity should not exceed the given order 4.
% This error means that in constructing your B-spline basis functions 
% (using a cubic spline, which has order 4), one of the knots (or input 
% data points) is duplicated too many times. In B-spline theory, the maximum 
% number of times (multiplicity) a knot can appear is limited by the order 
% of the spline. For a cubic B-spline (order 4), the knot multiplicity must 
% not exceed 4. 
% If your input data or knot vector contains a particular value more than 
% 4 times, spcol will throw this error.
% A common cause is having duplicate or nearly-identical data points when 
% generating the knot vector or control points. You should review your inputs 
% to ensure there are no unintended duplicates that inflate the knot 
% multiplicity.

% Even if your knots are distinct, you need to ensure that:
% 	•	The knot vector is properly constructed (including repeated endpoints as needed for a cubic spline),
% 	•	The evaluation points fall within the active support of the spline functions,
% 	•	And the spacing of the knots allows for nonzero basis functions over the domain of interest.
