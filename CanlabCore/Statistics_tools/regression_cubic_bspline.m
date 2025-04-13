function [out, residuals] = regression_cubic_bspline(X, y, varargin)
% regression_cubic_bspline - Fit a cubic B-spline regression model to one or more outcome variables with column expansion.
%
% :Usage:
% ::
%
%    [output_struct, residuals] = regression_cubic_bspline(X, y, [optional inputs])
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
%          OR
%         A matrix of multiple outcome values, with an outcome in each
%         column. A separate regression will be run for each column.
%         NOTE that the same spline basis is used, which will NOT be
%         optimal or appropriate unless the values are all on the same
%         scale (e.g., Z-scored). A future version may scale outcomes
%         internally so the same spline basis is appropriate and then
%         reverse the scaling back to the original scale after fitting. But
%         the current version requires the user to consider this and adjust
%         the input accordingly.
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
%    **names**: [cell array of strings] 
%         Names of regressors in the final design matrix.
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
%         - **F_test_table**: [table] F-statistics for each original regressor, with columns:
%               Regessor, F_stat, df1, df2, and p_value.
%               Rows are outcome variables in y, if y has more than 1 column.
%         - **P_values**: [table] P-value table for each original regressor, with columns:
%               Regessor, F_stat, df1, df2, and p_value.
%               Rows are outcome variables in y, if y has more than 1 column.
%         - **knot_pts**: [structure] Knot points used for each expanded column of X.
%         - **basis_sets**: [structure] B-spline basis matrices for each expanded column of X.
%         - **expanded_regs**: [vector] Indices of regressors (columns of X) that were expanded.
%         - **mapping_matrix**: [logical matrix]
%               Mapping from original regressors to columns of the expanded design matrix
%               (before adding the intercept).
%         - **input_params**: [structure] All input parameters passed to the function.
%         - **names**: [cell array of strings] Names of regressors in the final design matrix.
%         - **R2**: [1xK double] R-squared for each column of y (overall model fit).
%         - **R2_adj**: [1xK double] Adjusted R-squared for each column of y.
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
% :More examples of working with output:
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

addRequired(p, 'y', @(z) validateattributes(z, {'numeric'}, {'2d'}));

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

% validateattributes(y, {'numeric'}, {'numel', n_rows});
% The version below works for Y as a matrix
validateattributes(y, {'numeric'}, {'2d'});
assert(size(y, 1) == n_rows, 'y must have exactly %d rows', n_rows);

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

% -------------------------------------------------------------------------
% Build design matrix (X)
% -------------------------------------------------------------------------

% Expand Selected Columns using Cubic B-Spline Basis
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
            title('Basis function (linear)')
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


% Always add an intercept as the last column
% If it already has one somewhere in the middle, remove it and make sure it
% is at the end. Adjust names too.
wh = intercept(X_expanded, 'which');

if ~isempty(wh)
    X_expanded(wh, :) = [];
    new_names(wh, :) = [];
end

X_expanded = intercept(X_expanded, 'end');
new_names{end + 1} = 'Intercept';

%%

% -------------------------------------------------------------------------
% Fit Linear Regression Model
% -------------------------------------------------------------------------

beta  = X_expanded \ y;
y_hat = X_expanded * beta;

residuals = y - y_hat;

% -------------------------------------------------------------------------
% F-tests for each original regressor
% -------------------------------------------------------------------------
% Compute F-test by comparing full model to reduced model (excluding regressor i)
% Do this matrix-wise, so we can handle multiple columns of Y (multiple
% outcomes) simultaneously

% DFE and RSS for full model
n_obs = n_rows;
df2 = n_obs - size(X_expanded, 2);

RSS_full    = sum(residuals .^ 2, 1);  % Residual sums of squares for all outcomes

% Reduced models and stats

[X_reduced, beta_reduced, df1, RSS_reduced, F_stat, F_reg_names, p_values] = deal(cell(1, n_cols));  % reduced models for each outcome

for i = 1 : n_cols

    % Reduced models
    X_reduced{i} = X_expanded;
    X_reduced{i}(:, orig_idx{i}) = [];   % "k_i" = obs x reduced predictors for original column i

    % Model DF
    df1{i} = numel(orig_idx{i});  % error df_reduced - df_full

    % Betas
    beta_reduced{i} = X_reduced{i} \ y;  % k_i x num_outcomes; i is original design matrix column that has been (potentially) expanded

    % Reduced residual sums of squares
    RSS_reduced{i} = sum((y - X_reduced{i} * beta_reduced{i}).^2);

    % F-stats and P values
    F_stat{i} = ((RSS_reduced{i} - RSS_full) ./ df1{i}) ./ (RSS_full ./ df2);

    F_reg_names{i} = names{i};

    p_values{i} = 1 - fcdf(F_stat{i}, df1{i}, df2);

end

F_stats = cat(1, F_stat{:})'; % F stats, outcomes x original regressors

P_values = cat(1, p_values{:})'; % all P values, outcomes x original regressors

% adjust for P-values very close to zero
P_values(P_values < 100 * eps) = 100 * eps;

[partial_r2_mat, adj_r2_mat] = compute_partial_r2_sets(RSS_full, RSS_reduced, X_reduced, df2, n_rows);


% ---- Overall R² and Adjusted R² ----
y_mean = mean(y, 1);  % mean for each column
ss_total = sum((y - y_mean).^2, 1);  % total sum of squares per column

full_model_R2 = 1 - RSS_full ./ ss_total;
n = size(y, 1);
p = size(X, 2);  % includes intercept
R2_adj = 1 - (1 - full_model_R2) .* (n - 1) ./ (n - p);

% -------------------------------------------------------------------------
% Table for first outcome variable
% (typical use case is 1 outcome variable)
% -------------------------------------------------------------------------

df1s = cat(1, df1{:});
df2s = df2 * ones(size(df1s));

outcome_num = 1;

% for first outcome variable only
F_test_table = table(F_reg_names', F_stats(outcome_num, :)', df1s, df2s, P_values(outcome_num, :)', ...
    'VariableNames', {'Regressor', 'F_stat', 'df1', 'df2', 'p_value'});

% -------------------------------------------------------------------------
% Partial Regression Plots for Each Original Regressor and Actual vs.
% Predicted
% -------------------------------------------------------------------------
% First outcome only!
% Visualizes the unique contribution of one independent variable to the dependent variable, controlling for all other variables in the model.
% •	X-axis: Adjusted predictor: The residuals obtained from regressing the independent variable of interest X_i on all other independent variables.
% •	Y-axis: Adjusted data: The residuals obtained from regressing the dependent variable Y on all independent variables except X_i.


if input_params.doplot

    for i = 1 : n_cols

        % [y_res, X_i_res_orig, fit_line]
        plot_partial_residuals(X_expanded, X, y(:, outcome_num), orig_idx, i, names, fig_han);

    end % doplot

    % Plot Actual vs. Predicted if Requested

    create_figure('Design matrix and fits', 1, 2);

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

    scatter(y_hat(:, outcome_num), y(:, outcome_num), 'filled')
    lims = [min(min(y_hat(:, outcome_num)), min(y(:, outcome_num))), max(max(y_hat(:, outcome_num)), max(y(:, outcome_num)))];
    plot(lims, lims, 'r--', 'LineWidth', 1.5)

    hold off
    xlabel('Predicted y\_hat')
    ylabel('Actual y')
    title('Actual vs. Predicted (Cubic B-spline Regression)')

end  % plot


%% Assemble Output Structure


% Convert to tables
beta = array2table(beta', 'VariableNames', new_names);
F_stats = array2table(F_stats, 'VariableNames', input_params.names);
P_values = array2table(P_values, 'VariableNames', input_params.names);
partial_R2 = array2table(partial_r2_mat', 'VariableNames', input_params.names);
adj_R2 = array2table(adj_r2_mat', 'VariableNames', input_params.names);

out = struct();

out.X                 = X_expanded;
out.beta              = beta;
out.y_hat             = y_hat;
% out.partial_residuals = partial_residuals;

out.full_model_R2 = full_model_R2;
out.R2_adj = R2_adj;

out.F_stats = F_stats;
out.P_values = P_values;
out.partial_R2 = partial_R2;

out.adj_R2 = adj_R2;
out.F_test_table      = F_test_table;
out.knot_pts          = knot_pts;
out.basis_sets        = basis_sets;
out.expanded_regs     = expanded_regs;
out.mapping_matrix    = mapping_matrix;
out.was_expanded = sum(out.mapping_matrix, 2) > 1;
out.input_params      = input_params;
out.names             = new_names;


if input_params.verbose

    disp(' ')
    fprintf('Cubic B-spline regression completed. Model fitted with %d regressors.\n', size(X_expanded, 2));

    if size(y, 2) > 1  % multiple outcomes
        fprintf('Completed separate regressions for %3.0f outcomes (columns of y)\n', size(y, 2))
        disp('Showing regression table for first outcome: ')
        disp(' ')
    end

    disp(out.F_test_table)
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


function [partial_r2_mat, adj_r2_mat] = compute_partial_r2_sets(RSS_full, RSS_reduced, X_reduced, df2, n)
% Compute partial R² and adjusted partial R² across multiple outcomes and regressor sets.
%
% Inputs:
%   RSS_full        - 1 x k vector of residual sums of squares (full model), for k outcomes
%   RSS_reduced     - 1 x n cell, each a 1 x k vector of RSS for each reduced model
%   X_full          - [n x p_full] matrix: full model design matrix
%   X_reduced   - 1 x n cell, each a [n x p_i] matrix for reduced models
%   df1_all         - 1 x n vector of degrees of freedom added in each test
%   n               - number of observations
%
% Outputs:
%   partial_r2_mat  - [n x k] matrix of partial R² values
%   adj_r2_mat      - [n x k] matrix of adjusted partial R² values

% df1_all = cat(2, df1_all{:});     % Expand from cell array

n_sets = length(X_reduced);       % Number of regressor sets tested
k = length(RSS_full);             % Number of outcomes
% df2 = n - size(X_full, 2);        % df for full model

partial_r2_mat = zeros(n_sets, k);
adj_r2_mat = zeros(n_sets, k);

mse_full = RSS_full ./ df2;

for i = 1:n_sets

    rss_r = RSS_reduced{i};                 % 1 x k vector
    df0 = n - size(X_reduced{i}, 2);        % df for reduced model
    mse_reduced = rss_r / df0;

    partial_r2_mat(i, :) = (rss_r - RSS_full) ./ rss_r;
    adj_r2_mat(i, :)     = 1 - mse_full ./ mse_reduced;
end

end % subfunction


function [y_res, X_i_res_orig, fit_line] = plot_partial_residuals(X_expanded, X_original, y, orig_idx, i, names, fig_han)

% Subfunction to plot partial residual plots correctly

X_i = X_expanded(:, orig_idx{i});
X_i_orig = X_original(:, i);

X_other = X_expanded;
X_other(:, orig_idx{i}) = [];  % Remove the current predictor(s)

% Residuals of X_i regressed on other predictors
X_i_res = X_i - X_other * (X_other \ X_i);

% Use the non-expanded one for plot!
X_i_res_orig = X_i_orig - X_other * (X_other \ X_i_orig);

% Residuals of Y regressed on other predictors
y_res = y - X_other * (X_other \ y);

% Fit line: Regress Y residuals on X_i residuals
b_partial = X_i_res \ y_res;
fit_line = X_i_res * b_partial;

% Note: For a smooth estimate another option is to re-fit to partial residuals using the spline
% basis set.


% Plotting the partial residual plot
figure(fig_han(i))
subplot(1, 2, 2)
hold on
scatter(X_i_res_orig, y_res, 'filled')
[x_sorted, idx_sort] = sort(X_i_res_orig);
plot(x_sorted, fit_line(idx_sort), 'r-', 'LineWidth', 1.5)
hold off

xlabel([names{i} ' (Adjusted)'])
ylabel('Outcome (Adjusted)')
title(sprintf('Partial Regression Plot for %s', names{i}))

end
