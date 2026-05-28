function T = hrf_curve_summaries(source, varargin)
%HRF_CURVE_SUMMARIES Per-curve summary statistics from HRF score CSVs.
%
% Each row in a `*_map_scores.csv` is one (condition, lag) pair, and each
% score column is one signature, image-set map, or atlas region. The
% per-condition time course across lags IS the estimated HRF curve for
% that source. This helper extracts standard shape summaries (peak,
% time-to-peak, AUC, FWHM, onset, offset, duration, ...) per curve so
% curves can be directly compared across conditions, sources, subjects,
% runs, and models without re-deriving the same metrics ad hoc.
%
% Usage
% -----
%   T = hrf_curve_summaries(score_csv_path)
%   T = hrf_curve_summaries(score_table)
%   T = hrf_curve_summaries(input_table, 'TR', 0.8)
%
% Input dispatch (first arg)
% --------------------------
%   char/string ending in .csv  - load that single score CSV
%   table with 'volume_index'/'condition'/'lag_seconds' or 'lag_index'
%                               - treat as a single score table
%   table with 'subject', 'model', and beta_scores_file/t_scores_file
%                               - iterate rows of an input_table from
%                                 hrf_collect_wholebrain_outputs; load and
%                                 summarize every score CSV referenced
%
% Name-value parameters
% ---------------------
%   'Conditions'   - cellstr/string; subset to these condition names.
%                    Default {} (all conditions present in metadata).
%   'Sources'      - cellstr/string of column-name patterns (case-sensitive,
%                    partial-string match). Subsets which score columns
%                    are summarized. Default {} (all sig_*, map_*, atlas_*).
%   'Objects'      - which CSV objects to read from an input_table:
%                    {'beta','t'} or one of them. Default both.
%   'PeakThreshold'- fraction of |peak| used for onset/offset/FWHM.
%                    Default 0.5 (i.e., FWHM = full width at half max).
%   'TR'           - explicit TR if neither lag_seconds nor lag_index is
%                    populated. Default NaN.
%   'IncludeNaN'   - true to keep empty/all-NaN curves in the output (with
%                    metrics = NaN). Default false (skip them).
%
% Output
% ------
% T is a long table; one row per (file-of-origin, condition, source) curve:
%   subject, run_label, model, object       - origin metadata (when known)
%   source, source_kind, source_set, source_name
%                                            - 'sig_<set>_<name>' parsed
%   condition, n_lags, n_finite              - condition + curve coverage
%   peak_value, peak_lag_seconds, peak_lag_index
%   auc                                      - trapezoidal over lag_seconds
%   mean_amplitude, sd_amplitude
%   fwhm_seconds                             - full width at half peak
%   onset_lag_seconds, offset_lag_seconds, duration_seconds
%                                            - around the peak, at
%                                              PeakThreshold * |peak|
%
% Notes
% -----
% * peak is the *signed* value at argmax(|x|), so negative responses are
%   handled symmetrically with positive responses.
% * AUC is computed with trapz on lag_seconds; if all lags are equally
%   spaced this equals mean(x) * (max_lag - min_lag).
% * FWHM/onset/offset use the same |peak|-relative threshold so they're
%   self-consistent. For curves with a clean unimodal shape FWHM matches
%   the textbook definition; for noisy or multimodal curves FWHM is the
%   width of the *first contiguous* above-threshold run around the peak.
% * Score columns are detected by prefix: 'sig_*', 'map_*', 'atlas_*'.
%   SE / uncertainty columns (suffix '_se') are excluded.
%
% See also: hrf_score_one_prefix, hrf_score_wholebrain_input_table.

p = inputParser;
p.addRequired('source');
p.addParameter('Conditions', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Sources', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Objects', {'beta', 't'}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('PeakThreshold', 0.5, @(x) isscalar(x) && x > 0 && x < 1);
p.addParameter('TR', NaN, @(x) isscalar(x));
p.addParameter('IncludeNaN', false, @(x) islogical(x) || isnumeric(x));
p.parse(source, varargin{:});
opts = p.Results;

T = local_empty_summary_table();

if ischar(source) || isstring(source)
    score_table = local_read_csv(char(source));
    rows = local_summaries_one_table(score_table, opts);
    rows = local_attach_origin(rows, '', '', '', '');
    T = [T; rows];
    return
end

if istable(source)
    if local_looks_like_input_table(source)
        T = local_iterate_input_table(source, opts);
    else
        rows = local_summaries_one_table(source, opts);
        rows = local_attach_origin(rows, '', '', '', '');
        T = [T; rows];
    end
    return
end

error('hrf_curve_summaries:UnknownSource', ...
    'First argument must be a CSV path, a score table, or an input table. Got: %s.', class(source));
end


% =========================================================================
% Input-table dispatch
% =========================================================================
function tf = local_looks_like_input_table(t)
v = t.Properties.VariableNames;
has_object_paths = any(ismember(v, {'beta_scores_file', 't_scores_file', 'beta_file', 't_file'}));
tf = has_object_paths && any(ismember(v, {'subject', 'model'}));
end


function T = local_iterate_input_table(input_table, opts)
T = local_empty_summary_table();
object_kinds = lower(cellstr(string(opts.Objects)));

file_cols = struct('beta', 'beta_scores_file', 't', 't_scores_file');
for i = 1:height(input_table)
    subj = local_get_string(input_table, i, 'subject');
    run  = local_get_string(input_table, i, 'run_label');
    model = local_get_string(input_table, i, 'model');

    for k = 1:numel(object_kinds)
        obj = object_kinds{k};
        if ~isfield(file_cols, obj), continue; end
        col = file_cols.(obj);
        if ~any(strcmp(col, input_table.Properties.VariableNames)), continue; end
        path = char(string(input_table.(col)(i)));
        if isempty(path) || exist(path, 'file') ~= 2, continue; end
        try
            score_table = local_read_csv(path);
        catch err
            warning('hrf_curve_summaries:UnreadableCSV', ...
                'Skipping row %d %s: %s', i, path, err.message);
            continue
        end
        rows = local_summaries_one_table(score_table, opts);
        rows = local_attach_origin(rows, subj, run, model, obj);
        T = [T; rows]; %#ok<AGROW>
    end
end
end


% =========================================================================
% Per-table summarization
% =========================================================================
function T = local_summaries_one_table(score_table, opts)
T = local_empty_summary_table();
if isempty(score_table) || height(score_table) == 0
    return
end

[lag_seconds, lag_index, condition] = local_get_axes(score_table, opts.TR);
if isempty(lag_seconds) || isempty(condition)
    warning('hrf_curve_summaries:MissingAxes', ...
        'Score table is missing condition/lag axes; no curves summarized.');
    return
end

condition_list = local_filter_conditions(condition, opts.Conditions);
source_cols = local_score_columns(score_table, opts.Sources);

for c = 1:numel(condition_list)
    cond = condition_list{c};
    mask = strcmp(condition, cond);
    if ~any(mask), continue; end

    sub_lag_s = lag_seconds(mask);
    sub_lag_i = lag_index(mask);
    [sub_lag_s, ord] = sort(sub_lag_s);
    sub_lag_i = sub_lag_i(ord);

    for k = 1:numel(source_cols)
        col = source_cols{k};
        vals_raw = score_table.(col)(mask);
        vals = double(vals_raw(ord));
        metrics = local_compute_metrics(vals, sub_lag_s, opts.PeakThreshold);
        if ~opts.IncludeNaN && metrics.n_finite == 0
            continue
        end

        [skind, sset, sname] = local_parse_source(col);
        row = local_metric_row(col, skind, sset, sname, cond, sub_lag_i, sub_lag_s, metrics);
        T = [T; row]; %#ok<AGROW>
    end
end
end


function [lag_seconds, lag_index, condition] = local_get_axes(score_table, tr_arg)
v = score_table.Properties.VariableNames;
lag_seconds = [];
lag_index = [];
condition = strings(height(score_table), 1);

if any(strcmp('condition', v))
    condition = string(score_table.condition);
end

if any(strcmp('lag_seconds', v))
    lag_seconds = double(score_table.lag_seconds);
end
if any(strcmp('lag_index', v))
    lag_index = double(score_table.lag_index);
end

if isempty(lag_seconds) && ~isempty(lag_index) && ~isnan(tr_arg)
    lag_seconds = lag_index * tr_arg;
end
if isempty(lag_index) && ~isempty(lag_seconds)
    lag_index = lag_seconds;  % best we can do
end
if isempty(lag_seconds) && isempty(lag_index)
    return
end
end


function conds = local_filter_conditions(condition_vec, requested)
present = unique(cellstr(string(condition_vec)), 'stable');
if isempty(requested)
    conds = present(:)';
    return
end
requested = cellstr(string(requested));
conds = intersect(present, requested, 'stable');
conds = conds(:)';
end


function cols = local_score_columns(score_table, patterns)
v = score_table.Properties.VariableNames;
score_prefixes = {'sig_', 'map_', 'atlas_'};
keep = false(1, numel(v));
for i = 1:numel(v)
    name = v{i};
    if endsWith(name, '_se'), continue; end
    for p = 1:numel(score_prefixes)
        if startsWith(name, score_prefixes{p})
            keep(i) = true;
            break
        end
    end
end
cols = v(keep);

if ~isempty(patterns)
    patterns = cellstr(string(patterns));
    keep2 = false(size(cols));
    for i = 1:numel(cols)
        for j = 1:numel(patterns)
            if contains(cols{i}, patterns{j})
                keep2(i) = true;
                break
            end
        end
    end
    cols = cols(keep2);
end
end


% =========================================================================
% Source-name parsing  (sig_<set>_<name>, map_<set>_<name>, atlas_<set>_<name>_<suffix>)
% =========================================================================
function [kind, set_name, source_name] = local_parse_source(col)
kind = '';
set_name = '';
source_name = col;

parts = strsplit(col, '_');
if numel(parts) < 3, return; end

switch parts{1}
    case {'sig', 'map'}
        kind = parts{1};
        set_name = parts{2};
        source_name = strjoin(parts(3:end), '_');
    case 'atlas'
        kind = 'atlas';
        set_name = parts{2};
        % atlas columns end in a normalization suffix (_mean, _meanL1, _sum);
        % strip it from the source_name so the bare region label is visible.
        suffix_tokens = {'mean', 'meanL1', 'sum'};
        if ismember(parts{end}, suffix_tokens)
            source_name = strjoin(parts(3:end-1), '_');
        else
            source_name = strjoin(parts(3:end), '_');
        end
end
end


% =========================================================================
% Metric computation
% =========================================================================
function m = local_compute_metrics(vals, lag_seconds, peak_thresh_frac)
m = struct( ...
    'n_lags', numel(vals), ...
    'n_finite', sum(isfinite(vals)), ...
    'peak_value', NaN, ...
    'peak_lag_seconds', NaN, ...
    'peak_lag_index', NaN, ...
    'auc', NaN, ...
    'mean_amplitude', NaN, ...
    'sd_amplitude', NaN, ...
    'fwhm_seconds', NaN, ...
    'onset_lag_seconds', NaN, ...
    'offset_lag_seconds', NaN, ...
    'duration_seconds', NaN);

finite = isfinite(vals);
if ~any(finite)
    return
end

vfin = vals(finite);
lfin = lag_seconds(finite);

[~, k_peak_in_fin] = max(abs(vfin));
m.peak_value = vfin(k_peak_in_fin);
m.peak_lag_seconds = lfin(k_peak_in_fin);
% peak_lag_index = 1-based index into the original (unsorted-then-sorted) curve
peak_lag_pos = find(finite);
m.peak_lag_index = peak_lag_pos(k_peak_in_fin);

if abs(m.peak_value) > 0
    threshold = peak_thresh_frac * abs(m.peak_value);
    above = abs(vfin) >= threshold;
    % onset/offset are first/last lag in the *contiguous run containing the peak*
    [onset_k, offset_k] = local_contiguous_run(above, k_peak_in_fin);
    if ~isempty(onset_k)
        m.onset_lag_seconds = lfin(onset_k);
        m.offset_lag_seconds = lfin(offset_k);
        m.duration_seconds = lfin(offset_k) - lfin(onset_k);
        if peak_thresh_frac == 0.5
            m.fwhm_seconds = m.duration_seconds;
        else
            % If user passed a different threshold, still report FWHM at 0.5.
            half_above = abs(vfin) >= 0.5 * abs(m.peak_value);
            [hk0, hk1] = local_contiguous_run(half_above, k_peak_in_fin);
            if ~isempty(hk0)
                m.fwhm_seconds = lfin(hk1) - lfin(hk0);
            end
        end
    end
end

if numel(lfin) >= 2
    [lsort, ord] = sort(lfin);
    m.auc = trapz(lsort, vfin(ord));
end
m.mean_amplitude = mean(vfin);
m.sd_amplitude = std(vfin);
end


function [k0, k1] = local_contiguous_run(above, k_peak)
% Returns the [first, last] index of the contiguous run of true values in
% `above` that contains k_peak. If k_peak is not above-threshold, returns [].
k0 = [];
k1 = [];
if ~above(k_peak), return; end
% scan left
k0 = k_peak;
while k0 > 1 && above(k0 - 1)
    k0 = k0 - 1;
end
% scan right
k1 = k_peak;
while k1 < numel(above) && above(k1 + 1)
    k1 = k1 + 1;
end
end


% =========================================================================
% Row construction
% =========================================================================
function T = local_empty_summary_table()
T = table('Size', [0 18], ...
    'VariableTypes', { ...
        'string','string','string','string', ...
        'string','string','string','string', ...
        'string','double','double', ...
        'double','double','double', ...
        'double','double','double','double'}, ...
    'VariableNames', { ...
        'subject','run_label','model','object', ...
        'source','source_kind','source_set','source_name', ...
        'condition','n_lags','n_finite', ...
        'peak_value','peak_lag_seconds','peak_lag_index', ...
        'auc','mean_amplitude','sd_amplitude','fwhm_seconds'});
% extend with onset/offset/duration as separate vars (kept distinct for
% downstream filtering even though they live alongside fwhm)
T.onset_lag_seconds = double.empty(0, 1);
T.offset_lag_seconds = double.empty(0, 1);
T.duration_seconds = double.empty(0, 1);
end


function row = local_metric_row(col, skind, sset, sname, cond, ~, ~, metrics)
row = table( ...
    string(""), string(""), string(""), string(""), ...
    string(col), string(skind), string(sset), string(sname), ...
    string(cond), metrics.n_lags, metrics.n_finite, ...
    metrics.peak_value, metrics.peak_lag_seconds, metrics.peak_lag_index, ...
    metrics.auc, metrics.mean_amplitude, metrics.sd_amplitude, metrics.fwhm_seconds, ...
    'VariableNames', { ...
        'subject','run_label','model','object', ...
        'source','source_kind','source_set','source_name', ...
        'condition','n_lags','n_finite', ...
        'peak_value','peak_lag_seconds','peak_lag_index', ...
        'auc','mean_amplitude','sd_amplitude','fwhm_seconds'});
row.onset_lag_seconds = metrics.onset_lag_seconds;
row.offset_lag_seconds = metrics.offset_lag_seconds;
row.duration_seconds = metrics.duration_seconds;
end


function T = local_attach_origin(T, subject, run, model, object)
if height(T) == 0, return; end
T.subject(:) = string(subject);
T.run_label(:) = string(run);
T.model(:) = string(model);
T.object(:) = string(object);
end


% =========================================================================
% I/O helpers
% =========================================================================
function tbl = local_read_csv(path)
tbl = readtable(char(path), 'TextType', 'string');
end


function s = local_get_string(t, i, col)
s = '';
if any(strcmp(col, t.Properties.VariableNames))
    val = t.(col)(i);
    if isstring(val) || ischar(val)
        s = char(val);
    elseif iscell(val)
        s = char(val{1});
    else
        try
            s = char(string(val));
        catch
        end
    end
end
end
