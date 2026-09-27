function out = reliability(obj, varargin)
% reliability  ICC reliability of an rsm across replicates.
%
% Computes ICC (Shrout & Fleiss, 1979) on the vectorized upper triangle of
% the rsm across the replicate axis (3rd dim). The "items" are RSM cells
% (n_cells = k*(k-1)/2) and the "raters" are replicates (sessions, runs,
% or subjects, depending on R.level).
%
% Usage
% -----
%   % Single ICC across all replicates of the rsm
%   icc = R.reliability()
%   icc = R.reliability('icc_type','2-k')
%
%   % Per-subject ICC: compute ICC over each subject's replicates
%   T = R.reliability('by','subject_id')
%
% Optional name-value
% -------------------
%   'icc_type'  '3-k' (default) | '2-k' | '3-single' | '2-single' | '1-k' | '1-single'
%               Matches CanlabCore's ICC.m (case, type) syntax:
%                 '3-k' -> ICC(3,'k'), two-way mixed, average measures
%                 '2-k' -> ICC(2,'k'), two-way random, average measures
%               Average-measures variants use the average of the k replicates
%               and are typically what you want for RSM reliability across
%               replicates.
%
%   'by'        '' (default) or name of a replicate_table column for
%               per-group reliability (e.g. 'subject_id'). Returns a table
%               with one row per group.
%
%   'pool'      'auto' (default) | 'subject' | 'replicate'
%               How to aggregate ICC across the replicate axis when 'by' is
%               not set.
%                 'auto'     -> 'subject' if a subject column exists in
%                               replicate_table (subject_id / sub / etc.),
%                               else 'replicate'.
%                 'subject'  -> compute ICC PER SUBJECT across that subject's
%                               replicates (matches Sun et al. 01282025 paper
%                               recipe). Returns a struct with .summary
%                               (mean/median/std/n_subjects across subjects)
%                               and .per_subject (table).
%                 'replicate'-> pool all replicates as raters of the same
%                               items. Returns a scalar ICC. Conflates within
%                               and between subject variability -- use only
%                               when that's what you want.
%
%   'transform' 'auto' (default) | 'fisherz' | 'none'
%               'auto' applies atanh for correlation/spearman/cosine RSMs.
%
% Output
% ------
%   If 'by' is set:                table {group_var, ICC, n_replicates, n_cells}.
%   If pool='subject' (auto when subject col exists): struct with .summary
%                                  (mean/median/std/n_subjects) and .per_subject
%                                  table.
%   If pool='replicate':           scalar ICC.

if builtin('numel', obj) > 1
    error('rsm:reliability:nonScalar', 'expects a scalar rsm.');
end

p = inputParser;
p.addParameter('icc_type',           '3-k',  @(x) ischar(x) || isstring(x));
p.addParameter('by',                 '',     @(x) ischar(x) || isstring(x));
p.addParameter('transform',          'auto', @(x) ischar(x) || isstring(x));
p.addParameter('replicate_nan_frac', 0.5,    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('pool',               'auto', @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'auto','subject','replicate'}));
p.parse(varargin{:});
opt = p.Results;

% Resolve pool: 'auto' -> 'subject' if a subject column exists, else 'replicate'.
% Explicit 'by' overrides everything (it acts like per-group).
pool_mode = lower(char(opt.pool));
sub_col = '';
if isempty(opt.by) && ~strcmp(pool_mode, 'replicate')
    sub_col = detect_subject_column(obj.replicate_table);
    if isempty(sub_col) && strcmp(pool_mode, 'subject')
        warning('rsm:reliability:noSubjectColumn', ...
            ['pool=''subject'' requested but no subject column found in replicate_table. ', ...
             'Falling back to replicate-pooled ICC.']);
        pool_mode = 'replicate';
    elseif ~isempty(sub_col) && strcmp(pool_mode, 'auto')
        pool_mode = 'subject';
    else
        pool_mode = 'replicate';
    end
elseif ~isempty(opt.by)
    pool_mode = 'by';      % explicit per-group path -- skip auto logic
end

[case_id, type_id] = parse_icc_type(opt.icc_type);

% Vectorize each replicate's upper triangle
[k, ~, N] = size(obj.dat);
if N < 2
    error('rsm:reliability:tooFewReplicates', ...
        'Need at least 2 replicates (size(R,3) >= 2); got %d.', N);
end

upper = triu(true(k), 1);
n_cells = nnz(upper);
X = zeros(N, n_cells);
for n = 1:N
    slice = obj.dat(:, :, n);
    X(n, :) = slice(upper)';
end

X = apply_transform(X, opt.transform, obj);

if strcmp(pool_mode, 'replicate')
    % Pool across all replicates (raters = all replicates)
    out = compute_icc(X, case_id, type_id, opt.replicate_nan_frac);
    return
end

if strcmp(pool_mode, 'subject')
    % Per-subject ICC, then aggregate across subjects
    [G, subject_ids] = findgroups(obj.replicate_table.(sub_col));
    n_subj = max(G);
    icc_per_subj = nan(n_subj, 1);
    n_reps_per_subj = zeros(n_subj, 1);
    for g = 1:n_subj
        rows = find(G == g);
        n_reps_per_subj(g) = numel(rows);
        if numel(rows) < 2, continue; end
        icc_per_subj(g) = compute_icc(X(rows, :), case_id, type_id, opt.replicate_nan_frac);
    end
    % Return a struct with both per-subject and summary stats
    out = struct();
    out.pool = 'subject';
    out.summary = struct(...
        'mean',   mean(icc_per_subj, 'omitnan'), ...
        'median', median(icc_per_subj, 'omitnan'), ...
        'std',    std(icc_per_subj, 'omitnan'), ...
        'n_subjects', sum(~isnan(icc_per_subj)));
    out.per_subject = table(subject_ids, icc_per_subj, n_reps_per_subj, ...
        repmat(n_cells, n_subj, 1), ...
        'VariableNames', {char(sub_col), 'ICC', 'n_replicates', 'n_cells'});
    return
end

% Per-group ICC
by_col = char(opt.by);
if isempty(obj.replicate_table) || ~ismember(by_col, obj.replicate_table.Properties.VariableNames)
    error('rsm:reliability:noGroupCol', ...
        'replicate_table does not contain "%s".', by_col);
end

groups = obj.replicate_table.(by_col);
[G, group_ids] = findgroups(groups);
n_groups = max(G);
icc_vals = nan(n_groups, 1);
n_reps   = zeros(n_groups, 1);
for g = 1:n_groups
    rows = find(G == g);
    n_reps(g) = numel(rows);
    if numel(rows) < 2, continue; end
    icc_vals(g) = compute_icc(X(rows, :), case_id, type_id, opt.replicate_nan_frac);
end

out = table(group_ids, icc_vals, n_reps, repmat(n_cells, n_groups, 1), ...
    'VariableNames', {char(by_col), 'ICC', 'n_replicates', 'n_cells'});

end


function val = compute_icc(X, case_id, type_id, rep_nan_frac)
% X is [replicates x cells].
%
% Two-stage NaN handling:
%   1) Drop REPLICATES whose NaN fraction across cells exceeds rep_nan_frac
%      (default 0.5). These are typically "bad sessions" with whole missing
%      conditions producing NaN rows/cols in the per-replicate RSM.
%   2) Then drop CELLS that still have any NaN across the surviving replicates.
%
% This hybrid order recovers the across-all-replicates case where a minority
% of bad replicates would otherwise poison every cell. On clean data, both
% stages are no-ops and behavior matches the prior implementation.

if nargin < 4 || isempty(rep_nan_frac), rep_nan_frac = 0.5; end

% Stage 1: drop bad replicates
nan_frac_per_rep = mean(~isfinite(X), 2);
bad_reps = nan_frac_per_rep > rep_nan_frac;
if any(bad_reps)
    warning('rsm:reliability:droppedReplicates', ...
        'Dropped %d/%d replicates with > %.0f%% NaN cells before ICC.', ...
        sum(bad_reps), numel(bad_reps), 100 * rep_nan_frac);
    X = X(~bad_reps, :);
end
if size(X, 1) < 2
    warning('rsm:reliability:tooFewReplicates', ...
        ['After dropping bad replicates, only %d remain. Returning NaN. ', ...
         'Consider raising ''replicate_nan_frac'' or investigating data quality.'], ...
        size(X, 1));
    val = NaN; return
end

% Stage 2: drop NaN-containing cells across surviving replicates
data = X';   % cells x raters (ICC.m convention)
nan_rows = any(~isfinite(data), 2);
if any(nan_rows)
    n_drop = sum(nan_rows);
    n_keep = size(data, 1) - n_drop;
    if n_keep < 2
        warning('rsm:reliability:tooManyNaN', ...
            ['After dropping %d NaN-containing cells, only %d cells remain. ', ...
             'Returning NaN.'], n_drop, n_keep);
        val = NaN; return
    end
    data = data(~nan_rows, :);
end

if exist('ICC', 'file') == 2
    try
        val = ICC(case_id, type_id, data);
        if ~isnan(val), return; end
    catch
        % fall through to local
    end
end
val = local_icc(data, case_id, type_id);
end


function [case_id, type_id] = parse_icc_type(spec)
spec = lower(char(spec));
parts = strsplit(spec, '-');
if numel(parts) ~= 2
    error('rsm:reliability:badIccType', ...
        'icc_type must be e.g. ''3-k'', ''2-k'', ''3-single''; got %s.', spec);
end
case_id = str2double(parts{1});
type_id = parts{2};
if ~ismember(case_id, [1 2 3]) || ~ismember(type_id, {'k','single'})
    error('rsm:reliability:badIccType', ...
        'icc_type must be {1,2,3}-{single,k}; got %s.', spec);
end
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
switch mode
    case 'fisherz'
        X(X >  0.9999999) =  0.9999999;
        X(X < -0.9999999) = -0.9999999;
        X = atanh(X);
    case 'none'
        % no-op
end
end


function val = local_icc(data, case_id, type_id)
% Stock ICC fallback. Matches Shrout & Fleiss formulas for
% (case_id, type_id) in {1,2,3} x {single, k}.
[n, k] = size(data);    % n items, k raters
mean_items = mean(data, 2);
mean_raters = mean(data, 1);
grand = mean(data(:));

SS_t = sum((data(:) - grand).^2);
SS_b = k * sum((mean_items - grand).^2);             % between items
SS_w = SS_t - SS_b;                                   % within items
SS_r = n * sum((mean_raters - grand).^2);            % between raters
SS_e = SS_w - SS_r;                                   % residual (two-way)

MS_b = SS_b / (n - 1);
MS_w = SS_w / (n * (k - 1));
MS_r = SS_r / (k - 1);
MS_e = SS_e / ((n - 1) * (k - 1));

switch case_id
    case 1
        if strcmp(type_id, 'single'), val = (MS_b - MS_w) / (MS_b + (k - 1) * MS_w);
        else,                         val = (MS_b - MS_w) / MS_b;
        end
    case 2
        if strcmp(type_id, 'single')
            val = (MS_b - MS_e) / (MS_b + (k - 1) * MS_e + k * (MS_r - MS_e) / n);
        else
            val = (MS_b - MS_e) / (MS_b + (MS_r - MS_e) / n);
        end
    case 3
        if strcmp(type_id, 'single'), val = (MS_b - MS_e) / (MS_b + (k - 1) * MS_e);
        else,                         val = (MS_b - MS_e) / MS_b;
        end
end
end
