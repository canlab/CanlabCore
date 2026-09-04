function T = reliability_per_condition(obj, varargin)
% reliability_per_condition  Row-level ICC reliability — one per condition.
%
% For each condition (row of the RSM), pulls that row's relationships with
% all other conditions across replicates and computes ICC across replicates.
% Reproduces the per-individual-condition loop in `01282025 RSM Reliability`
% lines 144-204.
%
% Usage
% -----
%   T = R.reliability_per_condition()
%   T = R.reliability_per_condition('icc_type','2-k')
%
% Optional name-value
% -------------------
%   'icc_type'  '3-k' (default) | '2-k' | etc.
%   'by'        '' (default) | replicate_table column for per-(condition, by) ICC
%   'transform' 'auto' (default) | 'fisherz' | 'none'
%
% Output
% ------
%   T  table with one row per condition:
%        {Condition, ICC, n_replicates, n_other_conditions}

if builtin('numel', obj) > 1
    error('rsm:reliability_per_condition:nonScalar', 'expects a scalar rsm.');
end

p = inputParser;
p.addParameter('icc_type',           '3-k',  @(x) ischar(x) || isstring(x));
p.addParameter('by',                 '',     @(x) ischar(x) || isstring(x));
p.addParameter('transform',          'auto', @(x) ischar(x) || isstring(x));
p.addParameter('replicate_nan_frac', 0.5,    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('pool',               'auto', @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'auto','subject','replicate'}));
p.parse(varargin{:});
opt = p.Results;

% Resolve pool mode (matches reliability.m logic)
pool_mode = lower(char(opt.pool));
sub_col = '';
if isempty(opt.by) && ~strcmp(pool_mode, 'replicate')
    sub_col = detect_subject_column(obj.replicate_table);
    if isempty(sub_col) && strcmp(pool_mode, 'subject')
        warning('rsm:reliability_per_condition:noSubjectColumn', ...
            'pool=''subject'' requested but no subject column found; falling back to replicate.');
        pool_mode = 'replicate';
    elseif ~isempty(sub_col) && strcmp(pool_mode, 'auto')
        pool_mode = 'subject';
    else
        pool_mode = 'replicate';
    end
end

[case_id, type_id] = parse_icc_type(opt.icc_type);

[k, ~, N] = size(obj.dat);
if N < 2
    error('rsm:reliability_per_condition:tooFewReplicates', ...
        'Need at least 2 replicates; got %d.', N);
end

cond_labels = obj.labels;
if isempty(cond_labels)
    cond_labels = arrayfun(@(i) sprintf('cond_%d', i), (1:k)', 'UniformOutput', false);
end

% Pre-build per-condition X_c matrices
X_all = cell(k, 1);
for c = 1:k
    other_idx = [1:c-1, c+1:k];
    X_c = zeros(N, k - 1);
    for n = 1:N
        slice = obj.dat(:, :, n);
        X_c(n, :) = slice(c, other_idx);
    end
    X_all{c} = apply_transform(X_c, opt.transform, obj);
end

if strcmp(pool_mode, 'replicate')
    % Single ICC per condition pooled across all replicates
    icc_vals = nan(k, 1);
    for c = 1:k
        icc_vals(c) = compute_icc(X_all{c}, case_id, type_id, opt.replicate_nan_frac);
    end
    T = table(cond_labels, icc_vals, repmat(N, k, 1), repmat(k-1, k, 1), ...
        'VariableNames', {'Condition', 'ICC', 'n_replicates', 'n_other_conditions'});
    return
end

% pool='subject': per-condition × per-subject ICC, then aggregate
[G, subject_ids] = findgroups(obj.replicate_table.(sub_col));
n_subj = max(G);

mean_icc   = nan(k, 1);
median_icc = nan(k, 1);
std_icc    = nan(k, 1);
n_subj_v   = zeros(k, 1);
for c = 1:k
    icc_subj = nan(n_subj, 1);
    for g = 1:n_subj
        rows = find(G == g);
        if numel(rows) < 2, continue; end
        icc_subj(g) = compute_icc(X_all{c}(rows, :), case_id, type_id, opt.replicate_nan_frac);
    end
    mean_icc(c)   = mean(icc_subj,   'omitnan');
    median_icc(c) = median(icc_subj, 'omitnan');
    std_icc(c)    = std(icc_subj,    'omitnan');
    n_subj_v(c)   = sum(~isnan(icc_subj));
end

T = table(cond_labels, mean_icc, median_icc, std_icc, n_subj_v, ...
    repmat(k - 1, k, 1), ...
    'VariableNames', {'Condition', 'Mean_ICC', 'Median_ICC', 'Std_ICC', ...
                      'n_subjects', 'n_other_conditions'});

end


function val = compute_icc(X, case_id, type_id, rep_nan_frac)
% X is [replicates x cells]. See @rsm/reliability.m compute_icc for the
% two-stage NaN policy (drop bad replicates, then drop bad cells).
if nargin < 4 || isempty(rep_nan_frac), rep_nan_frac = 0.5; end
nan_frac_per_rep = mean(~isfinite(X), 2);
bad_reps = nan_frac_per_rep > rep_nan_frac;
if any(bad_reps)
    warning('rsm:reliability_per_condition:droppedReplicates', ...
        'Dropped %d/%d replicates with > %.0f%% NaN cells.', ...
        sum(bad_reps), numel(bad_reps), 100 * rep_nan_frac);
    X = X(~bad_reps, :);
end
if size(X, 1) < 2, val = NaN; return; end
data = X';
nan_rows = any(~isfinite(data), 2);
if any(nan_rows)
    if sum(~nan_rows) < 2, val = NaN; return; end
    data = data(~nan_rows, :);
end
if exist('ICC', 'file') == 2
    try, val = ICC(case_id, type_id, data); if ~isnan(val), return; end, catch, end
end
val = local_icc(data, case_id, type_id);
end


function [case_id, type_id] = parse_icc_type(spec)
spec = lower(char(spec));
parts = strsplit(spec, '-');
case_id = str2double(parts{1}); type_id = parts{2};
end


function X = apply_transform(X, mode, obj)
mode = lower(char(mode));
if strcmp(mode, 'auto')
    if obj.is_dissimilarity || ~ismember(lower(obj.metric), {'correlation','spearman','cosine'})
        mode = 'none';
    else
        mode = 'fisherz';
    end
end
if strcmp(mode, 'fisherz')
    X(X >  0.9999999) =  0.9999999;
    X(X < -0.9999999) = -0.9999999;
    X = atanh(X);
end
end


function val = local_icc(data, case_id, type_id)
[n, k] = size(data);
mean_items = mean(data, 2); mean_raters = mean(data, 1); grand = mean(data(:));
SS_t = sum((data(:) - grand).^2);
SS_b = k * sum((mean_items - grand).^2);
SS_w = SS_t - SS_b;
SS_r = n * sum((mean_raters - grand).^2);
SS_e = SS_w - SS_r;
MS_b = SS_b / (n - 1); MS_w = SS_w / (n * (k - 1));
MS_r = SS_r / (k - 1); MS_e = SS_e / ((n - 1) * (k - 1));
switch case_id
    case 1, if strcmp(type_id,'single'), val=(MS_b-MS_w)/(MS_b+(k-1)*MS_w); else, val=(MS_b-MS_w)/MS_b; end
    case 2, if strcmp(type_id,'single'), val=(MS_b-MS_e)/(MS_b+(k-1)*MS_e+k*(MS_r-MS_e)/n); else, val=(MS_b-MS_e)/(MS_b+(MS_r-MS_e)/n); end
    case 3, if strcmp(type_id,'single'), val=(MS_b-MS_e)/(MS_b+(k-1)*MS_e); else, val=(MS_b-MS_e)/MS_b; end
end
end
