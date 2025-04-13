function [X_windsorized, summary_table, was_windsorized] = windsorize_matrix_columnwise(X, varargin)
% Windsorize columns of a matrix or table to limit extreme values (outliers)
%
% :Usage:
% ::
%     [X_windsorized, summary_table, was_windsorized] = windsorize_matrix_columnwise(X, 'n_std_windsorize', 3, 'usemad', true);
%
% :Inputs:
%   **X:** [n x p] numeric matrix or table
%       Input data to Windsorize column-wise. If table, it must contain only numeric variables.
%
% :Optional Inputs:
%   **'n_std_windsorize', [numeric scalar]:**
%       Number of standard deviations (or MADs) to use as Windsorization threshold (default = 2.5).
%
%   **'doplot', [logical]:**
%       If true, generates histograms of original and Windsorized values (default = false).
%
%   **'verbose', [logical]:**
%       If true, prints summary table after processing (default = true).
%
%   **'usemad', [logical]:**
%       If true, use median absolute deviation (MAD) instead of standard deviation (default = false).
%
% :Outputs:
%   **X_windsorized:** matrix or table
%       Windsorized version of input X, returned in the same format.
%
%   **summary_table:** table
%       Summary statistics including mean, std/MAD, bounds, and percent Windsorized.
%
% :Examples:
% ::
%     T = table(randn(100,1)*10, randn(100,1)*5, 'VariableNames', {'A', 'B'});
%     [T_wind, summary] = windsorize_matrix_columnwise(T, 'n_std_windsorize', 2.5, 'usemad', true, 'verbose', true);
%
% :See also:
%   - zscore
%   - isoutlier
%

% Author: Your Name
% Copyright (C) 2024 Tor Wager

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

p = inputParser;
validate_input = @(x) isnumeric(x) && ismatrix(x) || istable(x);
p.addRequired('X', validate_input);
p.addParameter('n_std_windsorize', 3, @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0}));
p.addParameter('doplot', false, @(x) islogical(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.addParameter('usemad', false, @(x) islogical(x) && isscalar(x));
p.parse(X, varargin{:});
args = p.Results;

n_std = args.n_std_windsorize;
doplot = args.doplot;
verbose = args.verbose;
usemad = args.usemad;

% Convert table to matrix if needed
is_table_input = istable(X);
if is_table_input
    var_names = X.Properties.VariableNames;
    Xmat = table2array(X);
else
    Xmat = X;
    var_names = strcat("Var", string(1:size(Xmat, 2)));
end

[n, p] = size(Xmat);
X_windsorized = Xmat;
std_threshold = 1e-8;

% Initialize results
means = nan(p, 1);
scales = nan(p, 1);  % std or MAD
n_obs = nan(p, 1);
num_windsorized = nan(p, 1);
pct_windsorized = nan(p, 1);
lower_bounds = nan(p, 1);
upper_bounds = nan(p, 1);

if doplot
    n_rows = floor(sqrt(p));
    n_cols = ceil(p / n_rows);
    figure;
end

for i = 1:p
    col = Xmat(:, i);
    valid_idx = ~isnan(col);
    n_valid = sum(valid_idx);
    mu = mean(col(valid_idx));
    means(i) = mu;
    n_obs(i) = n_valid;

    if usemad
        scale = mad(col(valid_idx), 1);  % Use normalized MAD
    else
        scale = std(col(valid_idx));
    end
    scales(i) = scale;

    if scale < std_threshold
        lower_bounds(i) = NaN;
        upper_bounds(i) = NaN;
        pct_windsorized(i) = 0;
        continue;
    end

    lb = mu - n_std * scale;
    ub = mu + n_std * scale;
    lower_bounds(i) = lb;
    upper_bounds(i) = ub;

    col_w = col;
    windsor_mask = (col < lb | col > ub) & valid_idx;
    col_w(col < lb) = lb;
    col_w(col > ub) = ub;
    X_windsorized(:, i) = col_w;

    num_windsorized(i) = sum(windsor_mask);
    pct_windsorized(i) = 100 * num_windsorized(i) / n_valid;

    was_windsorized(:, i) = windsor_mask;

    if doplot
        subplot(n_rows, n_cols, i);

        nbins = ceil(size(X, 1) ./ 50);
        hh = histogram(col, nbins, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
        histogram(col_w, hh.BinEdges, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        title(var_names{i});
        if i == 1, legend({'Original', 'Windsorized'}); end
        axis tight;
    end

end % column



% Reformat result
if is_table_input
    X_windsorized = array2table(X_windsorized, 'VariableNames', var_names);
end

summary_table = table(var_names', means, scales, n_obs, num_windsorized, pct_windsorized, lower_bounds, upper_bounds, ...
    'VariableNames', {'variable', 'mean', 'mad_or_std', 'n_obs', 'num_windsorized', 'pct_windsorized', 'lower_bound', 'upper_bound'});

% Print summary table if verbose
if verbose
    disp(summary_table);
end

end