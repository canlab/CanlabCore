function G = hrf_curve_summary_groupstats(T, varargin)
%HRF_CURVE_SUMMARY_GROUPSTATS Group-level pooling of per-curve summaries.
%
% Takes the long table produced by hrf_curve_summaries and pools across a
% chosen axis (usually subjects), producing per-group descriptive
% statistics and a one-sample t-test of each metric against a null value.
%
% Usage
% -----
%   G = hrf_curve_summary_groupstats(T, ...
%           'GroupBy', {'condition','source_name'}, ...
%           'Across', 'subject', ...
%           'Metrics', {'peak_value','peak_lag_seconds','auc','fwhm_seconds'}, ...
%           'TestAgainst', 0, ...
%           'Correction', 'fdr');
%
% Inputs
% ------
%   T  -  table from hrf_curve_summaries (long form, one row per curve).
%
% Name-value parameters
% ---------------------
%   'GroupBy'      - cellstr of column names defining the groups (the
%                    "rows" of the output, one per (group, metric) cell).
%                    Default {'condition','source_name'}.
%   'Across'       - the column whose values are pooled within each group.
%                    Most commonly 'subject'. Default 'subject'. Set to ''
%                    to pool across ALL rows in each group.
%   'Metrics'      - cellstr of metric column names to test. Default = all
%                    numeric columns of T that aren't grouping/origin
%                    metadata.
%   'TestAgainst'  - null value for the one-sample t-test. Default 0.
%                    Pass NaN to skip the t-test (descriptives only).
%   'Alpha'        - significance threshold (used to populate the .sig
%                    boolean column and for FDR). Default 0.05.
%   'Correction'   - multiple-comparisons correction across rows of G:
%                      'none' (default), 'bonferroni', 'fdr'.
%   'MinN'         - minimum sample size to bother computing a t. Cells
%                    with n < MinN get NaN test stats. Default 3.
%
% Output
% ------
%   G is a long table with one row per (group cell x metric):
%       <GroupBy columns>, metric, n, mean, sd, sem, ci_lo, ci_hi,
%       null_value, t, df, p, p_corrected, sig
%   ci_lo/ci_hi are 95% confidence intervals from the t distribution.
%   sig is p_corrected < Alpha.
%
% Typical use
% -----------
%   % Per-subject curves -> group means and one-sample t per (condition,
%   % source_name) cell, FDR-corrected across rows:
%   G = hrf_curve_summary_groupstats(T_per_subject, ...
%           'GroupBy', {'condition','source_name'}, ...
%           'Across', 'subject', ...
%           'Correction', 'fdr');
%
%   % Subset to peak_lag for the NPS source, sorted by p:
%   Gnp = G(G.source_name == "NPS" & G.metric == "peak_lag_seconds", :);
%   sortrows(Gnp, 'p_corrected')
%
% Pairwise condition contrasts and unpaired group comparisons are NOT
% computed here; do those by pivoting T and running ttest() directly.
%
% See also: hrf_curve_summaries.

p = inputParser;
p.addRequired('T', @istable);
p.addParameter('GroupBy', {'condition', 'source_name'}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Across', 'subject', @(x) ischar(x) || isstring(x));
p.addParameter('Metrics', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('TestAgainst', 0, @(x) isscalar(x));
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.addParameter('Correction', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('MinN', 3, @(x) isscalar(x) && x >= 1);
p.parse(T, varargin{:});
opts = p.Results;

group_by = cellstr(string(opts.GroupBy));
across = char(opts.Across);
metrics = local_resolve_metrics(T, opts.Metrics);
null_value = opts.TestAgainst;
alpha = opts.Alpha;
correction = lower(strtrim(char(opts.Correction)));

G = local_empty_group_table(group_by, metrics);

if height(T) == 0 || isempty(metrics)
    return
end

% Build groupings.
group_keys = local_make_group_keys(T, group_by);
[unique_keys, ~, group_ix] = unique(group_keys, 'rows', 'stable');

for g = 1:size(unique_keys, 1)
    rows_g = T(group_ix == g, :);
    if ~isempty(across) && any(strcmp(across, T.Properties.VariableNames))
        % Collapse to one row per `across` value before pooling (defends
        % against the input having multiple runs per subject without prior
        % aggregation).
        rows_g = local_collapse_within_group(rows_g, across, metrics);
    end

    for m = 1:numel(metrics)
        metric = metrics{m};
        vals = local_get_metric_vals(rows_g, metric);
        row = local_groupstat_row(group_by, unique_keys(g, :), metric, vals, null_value, alpha, opts.MinN);
        G = [G; row]; %#ok<AGROW>
    end
end

% Across-row multiple-comparisons correction.
G = local_apply_correction(G, correction, alpha);
end


% =========================================================================
% Helpers
% =========================================================================
function metrics = local_resolve_metrics(T, requested)
if ~isempty(requested)
    metrics = cellstr(string(requested));
    metrics = metrics(:)';
    return
end
% Default: all numeric columns that aren't obvious origin/group fields.
exclude = {'subject','run_label','model','object','source','source_kind', ...
    'source_set','source_name','condition','n_lags','n_finite','n_modes', ...
    'n_sig_lags','peak_lag_index'};
metrics = {};
v = T.Properties.VariableNames;
for i = 1:numel(v)
    if ismember(v{i}, exclude), continue; end
    if isnumeric(T.(v{i}))
        metrics{end + 1} = v{i}; %#ok<AGROW>
    end
end
end


function keys = local_make_group_keys(T, group_by)
% Build an N x K string matrix of grouping-column values.
n = height(T);
k = numel(group_by);
keys = strings(n, k);
v = T.Properties.VariableNames;
for j = 1:k
    col = group_by{j};
    if ~any(strcmp(col, v))
        error('hrf_curve_summary_groupstats:UnknownGroupColumn', ...
            'GroupBy column "%s" is not in the input table.', col);
    end
    val = T.(col);
    if isstring(val)
        keys(:, j) = val;
    elseif iscell(val)
        keys(:, j) = string(val);
    else
        keys(:, j) = string(val);
    end
end
end


function rows_collapsed = local_collapse_within_group(rows_g, across, metrics)
% Take the mean of each metric per `across` value, so the t-test pools
% one value per subject (not one value per run/condition row).
if ~any(strcmp(across, rows_g.Properties.VariableNames))
    rows_collapsed = rows_g;
    return
end
ax_vals = rows_g.(across);
if isnumeric(ax_vals) || islogical(ax_vals)
    ax_vals = string(ax_vals);
elseif iscell(ax_vals)
    ax_vals = string(ax_vals);
end
[u, ~, ix] = unique(ax_vals, 'stable');
rows_collapsed = table();
rows_collapsed.(across) = u;
for m = 1:numel(metrics)
    if ~any(strcmp(metrics{m}, rows_g.Properties.VariableNames)), continue; end
    src = double(rows_g.(metrics{m}));
    out = NaN(numel(u), 1);
    for j = 1:numel(u)
        v = src(ix == j);
        out(j) = mean(v, 'omitnan');
    end
    rows_collapsed.(metrics{m}) = out;
end
end


function vals = local_get_metric_vals(rows_g, metric)
if ~any(strcmp(metric, rows_g.Properties.VariableNames))
    vals = [];
    return
end
vals = double(rows_g.(metric));
vals = vals(isfinite(vals));
end


function row = local_groupstat_row(group_by, key_vals, metric, vals, null_value, alpha, min_n)
n = numel(vals);
mu = NaN; sd = NaN; sem = NaN; ci_lo = NaN; ci_hi = NaN;
tval = NaN; df = NaN; pval = NaN; sig = false;
if n >= 1
    mu = mean(vals);
    if n >= 2
        sd = std(vals);
        sem = sd / sqrt(n);
    end
end
if isfinite(null_value) && n >= min_n && sd > 0 && isfinite(sd)
    df = n - 1;
    tval = (mu - null_value) / sem;
    pval = 2 * (1 - tcdf(abs(tval), df));
    tcrit = tinv(1 - alpha / 2, df);
    ci_lo = mu - tcrit * sem;
    ci_hi = mu + tcrit * sem;
    sig = pval < alpha;
end

row = table();
for j = 1:numel(group_by)
    row.(group_by{j}) = string(key_vals(j));
end
row.metric = string(metric);
row.n = n;
row.mean = mu;
row.sd = sd;
row.sem = sem;
row.ci_lo = ci_lo;
row.ci_hi = ci_hi;
row.null_value = null_value;
row.t = tval;
row.df = df;
row.p = pval;
row.p_corrected = pval;  % corrected later in a single pass
row.sig = sig;
end


function G = local_empty_group_table(group_by, metrics) %#ok<INUSD>
G = table();
for j = 1:numel(group_by)
    G.(group_by{j}) = strings(0, 1);
end
G.metric = strings(0, 1);
G.n = zeros(0, 1);
G.mean = zeros(0, 1);
G.sd = zeros(0, 1);
G.sem = zeros(0, 1);
G.ci_lo = zeros(0, 1);
G.ci_hi = zeros(0, 1);
G.null_value = zeros(0, 1);
G.t = zeros(0, 1);
G.df = zeros(0, 1);
G.p = zeros(0, 1);
G.p_corrected = zeros(0, 1);
G.sig = false(0, 1);
end


function G = local_apply_correction(G, correction, alpha)
if height(G) == 0, return; end
switch lower(correction)
    case {'none', '', 'unc', 'uncorrected'}
        % p_corrected = p (already initialized)
    case 'bonferroni'
        valid = isfinite(G.p);
        n_tests = sum(valid);
        if n_tests > 0
            G.p_corrected(valid) = min(G.p(valid) * n_tests, 1);
        end
    case 'fdr'
        valid = isfinite(G.p);
        if any(valid)
            p_in = G.p(valid);
            [p_sorted, ord] = sort(p_in);
            n_tests = numel(p_sorted);
            ranks = (1:n_tests)';
            bh = p_sorted .* (n_tests ./ ranks);
            for k = n_tests - 1 : -1 : 1
                bh(k) = min(bh(k), bh(k + 1));
            end
            bh = min(bh, 1);
            % Re-order back to original positions
            p_corr_back = NaN(n_tests, 1);
            p_corr_back(ord) = bh;
            G.p_corrected(valid) = p_corr_back;
        end
    otherwise
        warning('hrf_curve_summary_groupstats:UnknownCorrection', ...
            'Unknown Correction: %s. Using none.', correction);
end
G.sig = G.p_corrected < alpha;
end
