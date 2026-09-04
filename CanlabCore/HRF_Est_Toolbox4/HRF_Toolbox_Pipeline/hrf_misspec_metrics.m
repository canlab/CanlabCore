function T = hrf_misspec_metrics(source, varargin)
%HRF_MISSPEC_METRICS Curve-shape misspecification metrics vs a reference HRF.
%
% For each estimated HRF curve (one per subject/run/model/condition/source
% column in the score CSVs), quantifies how its SHAPE departs from a
% reference HRF -- the canonical SPM HRF by default, or an empirical /
% user-supplied reference. This is the curve-based half of the Phase 4
% misspecification pipeline; residual-based metrics (ResidScan, PowerLoss)
% live separately because they need the residual time course, which is not
% stored in the whole-brain score CSVs.
%
% Usage
% -----
%   T = hrf_misspec_metrics(input_table)                       % vs canonical
%   T = hrf_misspec_metrics(input_table, 'Source','signature','Set','all')
%   T = hrf_misspec_metrics(input_table, 'Reference','custom', ...
%           'ReferenceHRF', h_vec, 'ReferenceLags', t_vec)
%   T = hrf_misspec_metrics(struct('label',{'acc','exp'},'table',{it1,it2}))
%
% Input dispatch (first arg) -- same shapes as hrf_curve_summaries:
%   char/string ending .csv  -> one score CSV
%   table with condition+lag  -> a single score table
%   input_table (subject/model/*_scores_file) -> iterate its rows
%   struct array .label/.table -> per-source, rows tagged with study_label
%
% Name-value parameters
% ---------------------
%   'Reference'    - 'canonical' (default, SPM spm_hrf), 'empirical', or
%                    'custom'.
%   'ReferenceHRF' - for 'custom': vector of HRF values. With
%                    'ReferenceLags' it is resampled to each curve's lags;
%                    without, it must already match the curve length.
%   'ReferenceLags'- for 'custom': lag-seconds axis of ReferenceHRF.
%   'EmpiricalRef' - for 'empirical': a containers.Map or struct keyed by
%                    condition (or 'all') -> reference HRF vector, plus
%                    'EmpiricalRefLags'. (Empirical group-optimal HRF from
%                    HMHRFest is a v1 deliverable; supply your own here.)
%   'Source'/'Set' - column family selection. 'atlas' (default),
%                    'signature', 'imageset', a cellstr of these, or 'all'
%                    (atlas + signature + imageset in one table). Within
%                    'atlas', ALL atlases present are returned (source_set =
%                    the atlas token, e.g. canlab2024 vs painpathways2024).
%   'Set'          - which SET to keep, matched against source_set: the atlas
%                    token ('canlab2024'), signature set ('all'), or imageset
%                    ('bucknerlab'). char/string/cellstr + glob wildcards.
%                    Default {} = all sets.
%   'Names'        - which specific source to keep, matched against source_name:
%                    the region ('PAG','*Ins*'), signature ('NPS'), or network
%                    map name. char/string/cellstr + glob wildcards. Default
%                    {} = all. Set picks the atlas/set; Names picks the item.
%   'Conditions'   - cellstr; subset (glob wildcards ok). Default all.
%   'Objects'      - {'beta','t'} for input_table iteration. Default {'beta'}.
%   'Model'        - filter input_table rows by model. Default '' (all).
%   'TR'           - fallback if a CSV lacks lag_seconds. Default NaN.
%   'MinLags'      - skip curves with fewer than this many finite lags.
%                    Default 4.
%   'GroupCurveFirst' - false (default): one metric row per subject/run curve
%                    (then pool with hrf_curve_summary_groupstats). true:
%                    average the curves ACROSS subjects/runs within each
%                    (model, object) FIRST, then score the group-mean curve.
%                    Far less noise-sensitive -- use it when per-subject HRFs
%                    are noisy (sparse conditions, low SNR). Output rows carry
%                    subject='group' and need NO groupstats.
%
% Output
% ------
% Long table; one row per (origin, condition, source) curve:
%   subject, run_label, model, object, study_label (when multi-source)
%   source, source_kind, source_set, source_name, condition, n_lags
%   reference                         - which reference was used
%   peak_lag_seconds, ref_peak_lag_seconds, peak_lag_bias_seconds
%   auc, ref_auc, auc_ratio
%   shape_corr                        - Pearson corr(curve, ref) over lags
%   shape_r2                          - shape_corr^2 (variance shared with ref)
%   misspec_r2                        - 1 - SSresid/SStot after best-scaling
%                                       the reference to the curve; LOW or
%                                       negative => poorly described by the
%                                       reference shape => more misspecified
%   scaled_resid_rmse                 - rmse(curve - a*ref), a = LS scale
%
% Notes
% -----
% * shape_corr / shape_r2 are scale- and sign-invariant measures of how
%   well the curve tracks the reference time course. misspec_r2 additionally
%   accounts for amplitude via a least-squares scale of the reference.
% * peak_lag_bias > 0 means the estimated HRF peaks LATER than the reference.
% * auc_ratio > 1 means more integrated response than the reference.
%
% See also: hrf_curve_summaries, plot_hrf_curves, ResidScan, PowerLoss.

p = inputParser;
p.addRequired('input_source');   % label must differ from the 'Source' param
p.addParameter('Reference', 'canonical', @(x) ischar(x) || isstring(x));
p.addParameter('ReferenceHRF', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('ReferenceLags', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('EmpiricalRef', [], @(x) isempty(x) || isstruct(x) || isa(x, 'containers.Map'));
p.addParameter('EmpiricalRefLags', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('Source', 'atlas', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('Set', '', @(x) ischar(x) || isstring(x) || iscell(x));
% Names: which specific source(s) to keep, matched against the parsed source
% NAME (region / signature / imageset-map name), cellstr/string + glob
% wildcards. Default {} = all. (Set selects the atlas/signature-set/imageset;
% Names selects the item within it.)
p.addParameter('Names', {}, @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('Conditions', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Objects', {'beta'}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Model', '', @(x) ischar(x) || isstring(x));
p.addParameter('TR', NaN, @(x) isscalar(x));
p.addParameter('MinLags', 4, @(x) isscalar(x) && x >= 2);
% GroupCurveFirst: average the HRF curves ACROSS subjects/runs (within each
% model+object) before computing the misspecification metrics, instead of
% computing per-subject metrics that a later groupstats averages. The group-
% mean curve is far less noise-sensitive, so misspec_r2 / peak_lag_bias
% reflect the true group HRF shape rather than per-subject estimation noise.
% Output rows then carry subject='group' and need NO hrf_curve_summary_groupstats.
p.addParameter('GroupCurveFirst', false, @(x) islogical(x) || isnumeric(x));
p.parse(source, varargin{:});
opts = p.Results;

T = local_empty_metrics_table();

% Multi-source struct array: recurse, tag study_label.
if isstruct(source) && isfield(source, 'label') && ...
        (isfield(source, 'table') || isfield(source, 'input_table'))
    chunks = cell(numel(source), 1);
    for i = 1:numel(source)
        if isfield(source(i), 'table'), sub = source(i).table; else, sub = source(i).input_table; end
        Ti = hrf_misspec_metrics(sub, varargin{:});
        if ~isempty(Ti) && height(Ti) > 0
            Ti.study_label = repmat(string(source(i).label), height(Ti), 1);
        end
        chunks{i} = Ti;
    end
    chunks = chunks(~cellfun(@isempty, chunks));
    if ~isempty(chunks), T = vertcat(chunks{:}); end
    return
end

% Resolve the table list to process.
if ischar(source) || isstring(source)
    score_tables = {readtable(char(source), 'TextType', 'string')};
    origins = struct('subject', {''}, 'run_label', {''}, 'model', {''}, 'object', {''});
elseif istable(source) && local_looks_like_input_table(source)
    [score_tables, origins] = local_load_input_table(source, opts);
elseif istable(source)
    score_tables = {source};
    origins = struct('subject', {''}, 'run_label', {''}, 'model', {''}, 'object', {''});
else
    error('hrf_misspec_metrics:UnknownSource', ...
        'First arg must be a CSV path, a score table, or an input_table.');
end

% Source can be one family ('atlas'/'signature'/'imageset'), a cellstr of
% families, or 'all' (= all three). Each family is matched separately and the
% rows concatenated; source_kind/source_set in the output distinguish them
% (and, for atlas, source_set is the atlas token so canlab2024 vs
% painpathways2024 are kept apart).
families = local_resolve_families(opts.Source);

% GroupCurveFirst: replace the per-subject tables with one group-mean curve
% table per (model, object), so metrics are computed on the averaged HRF.
if logical(opts.GroupCurveFirst) && numel(score_tables) > 1
    [score_tables, origins] = local_group_curve_pool(score_tables, origins);
end

rows = {};
for ti = 1:numel(score_tables)
    St = score_tables{ti};
    if isempty(St) || height(St) == 0, continue; end
    org = origins(ti);
    for f = 1:numel(families)
        ms = struct('family', families{f}, 'set', char(opts.Set));
        chunk = local_metrics_one_table(St, org, ms, opts);
        if ~isempty(chunk), rows{end + 1} = chunk; end %#ok<AGROW>
    end
end
if ~isempty(rows)
    T = vertcat(rows{:});
end
end


% =========================================================================
% Per-table metric computation
% =========================================================================
function T = local_metrics_one_table(St, org, match_spec, opts)
T = local_empty_metrics_table();
v = St.Properties.VariableNames;
if ~any(strcmp('condition', v)) || ~any(strcmp('lag_seconds', v))
    if ~any(strcmp('lag_seconds', v)) && any(strcmp('lag_index', v)) && ~isnan(opts.TR)
        St.lag_seconds = double(St.lag_index) * opts.TR;
    else
        return
    end
end

cond_all = string(St.condition);
cond_list = local_filter_conditions(cond_all, opts.Conditions);
[score_cols, score_labels] = local_match_columns(v, match_spec);
if isempty(score_cols), return; end

rows = {};
for c = 1:numel(cond_list)
    cmask = cond_all == string(cond_list{c});
    if ~any(cmask), continue; end
    lag = double(St.lag_seconds(cmask));
    [lag, ord] = sort(lag);

    for k = 1:numel(score_cols)
        [skind, sset, sname] = local_parse_source(score_cols{k}, score_labels{k}, match_spec);
        if ~local_match_any(sset, opts.Set), continue; end       % which atlas / sig-set / imageset
        if ~local_match_any(sname, opts.Names), continue; end    % which specific region / signature / map
        y = double(St.(score_cols{k})(cmask));
        y = y(ord);
        fin = isfinite(y) & isfinite(lag);
        if sum(fin) < opts.MinLags, continue; end
        yy = y(fin); xx = lag(fin);

        ref = local_reference_hrf(xx, cond_list{c}, opts);
        if isempty(ref) || all(~isfinite(ref)), continue; end

        m = local_compute_misspec(yy, ref, xx);
        rows{end + 1} = local_metric_row(org, score_cols{k}, skind, sset, sname, ...
            cond_list{c}, char(opts.Reference), m); %#ok<AGROW>
    end
end
if ~isempty(rows)
    T = vertcat(rows{:});
end
end


function [pooled_tables, pooled_origins] = local_group_curve_pool(score_tables, origins)
% Group the per-subject score tables by (model, object) and average each
% group's curves into one group-mean table -- so GroupCurveFirst computes the
% metrics on the averaged HRF rather than per subject.
models  = local_origin_field(origins, 'model');
objects = local_origin_field(origins, 'object');
key = strcat(models, "||", objects);
[ukey, ~, gidx] = unique(key, 'stable');
pooled_tables = cell(1, numel(ukey));
pooled_origins = repmat(struct('subject', 'group', 'run_label', '', 'model', '', 'object', ''), 1, numel(ukey));
for u = 1:numel(ukey)
    idx = find(gidx == u);
    pooled_tables{u} = local_pool_curves_across_tables(score_tables(idx));
    o = origins(idx(1));
    pooled_origins(u) = struct('subject', 'group', 'run_label', '', ...
        'model', local_origin_value(o, 'model'), 'object', local_origin_value(o, 'object'));
end
end

function s = local_origin_field(origins, f)
s = strings(1, numel(origins));
for i = 1:numel(origins)
    s(i) = string(local_origin_value(origins(i), f));
end
end

function v = local_origin_value(o, f)
if isfield(o, f) && ~isempty(o.(f)), v = char(string(o.(f))); else, v = ''; end
end

function P = local_pool_curves_across_tables(tables)
% Average numeric columns within (condition, lag) across same-structured score
% tables -> one group-mean curve per condition x lag.
%
% Robust to PARTIAL tables: instead of keeping only the column intersection
% (which lets a single stale/failed score CSV that is missing the atlas/sig/
% map columns silently collapse the pool to whatever few columns every table
% happens to share), this takes the UNION of numeric columns and averages each
% one over only the tables that actually contain it (omitnan). A subject
% missing canlab2024 thus drops out of canlab2024's group mean rather than
% deleting canlab2024 for everyone. Uneven coverage raises a warning so the
% offending CSVs get re-scored.
tables = tables(~cellfun(@isempty, tables));
if isempty(tables), P = table(); return; end
if isscalar(tables), P = tables{1}; return; end

% Group key: condition + a lag key present in EVERY table.
have_lagidx = all(cellfun(@(t) any(strcmp('lag_index', t.Properties.VariableNames)), tables));
if have_lagidx, lagkey = 'lag_index'; else, lagkey = 'lag_seconds'; end
for i = 1:numel(tables)
    vn = tables{i}.Properties.VariableNames;
    if ~any(strcmp('condition', vn)) || ~any(strcmp(lagkey, vn))
        error('hrf_misspec_metrics:PoolKeys', ...
            'GroupCurveFirst pooling needs ''condition'' and ''%s'' in every score table.', lagkey);
    end
end

% Union of numeric value columns (everything except the two group keys),
% in first-seen order.
valnames = {};
for i = 1:numel(tables)
    ti = tables{i};
    vn = ti.Properties.VariableNames;
    for j = 1:numel(vn)
        nm = vn{j};
        if strcmp(nm, 'condition') || strcmp(nm, lagkey), continue; end
        if ~isnumeric(ti.(nm)), continue; end
        if ~any(strcmp(nm, valnames)), valnames{end + 1} = nm; end %#ok<AGROW>
    end
end

% Stack, NaN-filling columns a given table lacks; track per-column coverage.
parts = cell(1, numel(tables));
covcount = zeros(1, numel(valnames));
for i = 1:numel(tables)
    ti = tables{i};
    h = height(ti);
    S = table();
    S.condition = string(ti.condition);
    S.(lagkey) = double(ti.(lagkey));
    has = ismember(valnames, ti.Properties.VariableNames);
    covcount = covcount + has;
    for j = 1:numel(valnames)
        if has(j) && isnumeric(ti.(valnames{j}))
            S.(valnames{j}) = double(ti.(valnames{j}));
        else
            S.(valnames{j}) = nan(h, 1);
        end
    end
    parts{i} = S;
end
big = vertcat(parts{:});

[g, gcond, glag] = findgroups(big.condition, big.(lagkey));
P = table();
P.condition = gcond;
P.(lagkey) = glag;
for j = 1:numel(valnames)
    P.(valnames{j}) = splitapply(@(x) mean(x, 'omitnan'), big.(valnames{j}), g);
end

under = covcount < numel(tables);
if any(under)
    warning('hrf_misspec_metrics:PartialCoverage', ...
        ['GroupCurveFirst: %d of %d pooled columns are absent from some score tables ' ...
         '(as few as %d/%d tables carry a column); each is averaged over only the tables ' ...
         'that have it. This usually means stale/partial score CSVs -- run ' ...
         'hrf_audit_score_freshness on the output dir and re-score the flagged runs. ' ...
         '%d columns are fully covered.'], ...
        sum(under), numel(valnames), min(covcount(under)), numel(tables), sum(~under));
end
end


function m = local_compute_misspec(y, ref, x)
m = struct('n_lags', numel(y), ...
    'peak_lag_seconds', NaN, 'ref_peak_lag_seconds', NaN, 'peak_lag_bias_seconds', NaN, ...
    'auc', NaN, 'ref_auc', NaN, 'auc_ratio', NaN, ...
    'shape_corr', NaN, 'shape_r2', NaN, 'misspec_r2', NaN, 'scaled_resid_rmse', NaN);

% Peak lags (signed argmax of magnitude).
[~, ky] = max(abs(y));   m.peak_lag_seconds = x(ky);
[~, kr] = max(abs(ref)); m.ref_peak_lag_seconds = x(kr);
m.peak_lag_bias_seconds = m.peak_lag_seconds - m.ref_peak_lag_seconds;

% AUC.
if numel(x) >= 2
    m.auc = trapz(x, y);
    m.ref_auc = trapz(x, ref);
    if m.ref_auc ~= 0, m.auc_ratio = m.auc / m.ref_auc; end
end

% Shape correlation.
if std(y) > 0 && std(ref) > 0
    r = corr(y(:), ref(:));
    m.shape_corr = r;
    m.shape_r2 = r ^ 2;
end

% Best least-squares scale of the reference, then residual R2 + rmse.
denom = ref(:)' * ref(:);
if denom > 0
    a = (ref(:)' * y(:)) / denom;
    resid = y(:) - a * ref(:);
    ss_tot = sum((y(:) - mean(y)) .^ 2);
    if ss_tot > 0
        m.misspec_r2 = 1 - sum(resid .^ 2) / ss_tot;
    end
    m.scaled_resid_rmse = sqrt(mean(resid .^ 2));
end
end


% =========================================================================
% Reference HRF resolution
% =========================================================================
function ref = local_reference_hrf(lags, condition, opts)
ref = [];
switch lower(strtrim(char(opts.Reference)))
    case 'canonical'
        ref = local_canonical_at_lags(lags);
    case 'custom'
        if isempty(opts.ReferenceHRF)
            error('hrf_misspec_metrics:NoCustomRef', ...
                'Reference=''custom'' requires ''ReferenceHRF''.');
        end
        ref = local_resample_ref(opts.ReferenceHRF, opts.ReferenceLags, lags);
    case 'empirical'
        if isempty(opts.EmpiricalRef)
            error('hrf_misspec_metrics:NoEmpiricalRef', ...
                ['Reference=''empirical'' requires ''EmpiricalRef'' (per-condition ' ...
                 'HRF vectors). The group-optimal empirical HRF from HMHRFest is a ' ...
                 'v1 deliverable -- for now pass your own via EmpiricalRef + ' ...
                 'EmpiricalRefLags.']);
        end
        rv = local_empirical_lookup(opts.EmpiricalRef, condition);
        ref = local_resample_ref(rv, opts.EmpiricalRefLags, lags);
    otherwise
        error('hrf_misspec_metrics:UnknownReference', ...
            'Unknown Reference: %s. Use canonical, empirical, or custom.', char(opts.Reference));
end
end


function ref = local_canonical_at_lags(lags)
ref = [];
if exist('spm_hrf', 'file') ~= 2
    error('hrf_misspec_metrics:NoSPM', ...
        'Reference=''canonical'' needs spm_hrf on the path (SPM).');
end
dt = 0.1;
h = spm_hrf(dt);
h = h ./ max(h);
t_fine = (0:numel(h) - 1) * dt;
ref = interp1(t_fine, h, lags(:), 'linear', 0);
end


function out = local_resample_ref(vals, ref_lags, lags)
vals = double(vals(:));
if isempty(ref_lags)
    if numel(vals) == numel(lags)
        out = vals;
    else
        error('hrf_misspec_metrics:RefLengthMismatch', ...
            ['Reference vector length (%d) does not match the curve length (%d) ' ...
             'and no ReferenceLags were given to resample.'], numel(vals), numel(lags));
    end
    return
end
ref_lags = double(ref_lags(:));
out = interp1(ref_lags, vals, lags(:), 'linear', 0);
end


function rv = local_empirical_lookup(empref, condition)
rv = [];
if isa(empref, 'containers.Map')
    if isKey(empref, char(condition)), rv = empref(char(condition));
    elseif isKey(empref, 'all'), rv = empref('all'); end
elseif isstruct(empref)
    f = matlab.lang.makeValidName(char(condition));
    if isfield(empref, f), rv = empref.(f);
    elseif isfield(empref, 'all'), rv = empref.all; end
end
if isempty(rv)
    error('hrf_misspec_metrics:EmpiricalRefMissing', ...
        'EmpiricalRef has no entry for condition ''%s'' or ''all''.', char(condition));
end
end


% =========================================================================
% Column matching (mirrors plot_hrf_curves families)
% =========================================================================
function fams = local_resolve_families(src)
% Expand the Source argument into a cellstr of column families.
s = lower(strtrim(string(src)));
s = s(strlength(s) > 0);
if isempty(s), s = "atlas"; end
if any(s == "all"), fams = {'atlas', 'signature', 'imageset'}; return; end
fams = cellstr(s);
end

function [cols, labels] = local_match_columns(v, match_spec)
cols = {}; labels = {};
switch match_spec.family
    case 'atlas',     prefix = 'atlas_';
    case 'signature', prefix = 'sig_';
    case 'imageset',  prefix = 'map_';
    otherwise,        prefix = 'atlas_';
end
% Match every column of the family; Set/Names selection is applied later on
% the PARSED source identity (source_set / source_name) in the caller.
for i = 1:numel(v)
    name = v{i};
    if ~startsWith(name, prefix), continue; end
    if endsWith(name, '_se'), continue; end
    cols{end + 1} = name; %#ok<AGROW>
    labels{end + 1} = name; %#ok<AGROW>
end
end

function tf = local_match_any(value, patterns)
% True if VALUE matches any of PATTERNS (glob wildcards * ?; exact otherwise).
% Empty PATTERNS => match-all (no filter).
patterns = cellstr(string(patterns));
patterns = patterns(~cellfun(@(s) isempty(strtrim(s)), patterns));
if isempty(patterns), tf = true; return; end
val = char(string(value));
tf = false;
for p = 1:numel(patterns)
    pat = strtrim(patterns{p});
    if any(pat == '*' | pat == '?')
        if ~isempty(regexp(val, ['^', regexptranslate('wildcard', pat), '$'], 'once'))
            tf = true; return
        end
    elseif strcmpi(val, pat)
        tf = true; return
    end
end
end


function [kind, set_name, source_name] = local_parse_source(col, ~, match_spec)
kind = match_spec.family;
set_name = '';
source_name = col;
parts = strsplit(col, '_');
if numel(parts) >= 3
    set_name = parts{2};
    switch match_spec.family
        case {'signature', 'imageset'}
            source_name = strjoin(parts(3:end), '_');
        case 'atlas'
            suffix_tokens = {'mean', 'meanL1', 'sum'};
            if ismember(parts{end}, suffix_tokens)
                source_name = strjoin(parts(3:end-1), '_');
            else
                source_name = strjoin(parts(3:end), '_');
            end
    end
end
end


% =========================================================================
% Input-table handling
% =========================================================================
function tf = local_looks_like_input_table(t)
v = t.Properties.VariableNames;
tf = any(ismember(v, {'beta_scores_file', 't_scores_file'})) && ...
    any(ismember(v, {'subject', 'model'}));
end


function [tables, origins] = local_load_input_table(input_table, opts)
tables = {};
origins = struct('subject', {}, 'run_label', {}, 'model', {}, 'object', {});
object_kinds = lower(cellstr(string(opts.Objects)));
file_cols = struct('beta', 'beta_scores_file', 't', 't_scores_file');
model_filter = lower(char(opts.Model));
for i = 1:height(input_table)
    row_model = local_get_string(input_table, i, 'model');
    if ~isempty(model_filter) && ~strcmpi(row_model, model_filter), continue; end
    for k = 1:numel(object_kinds)
        obj = object_kinds{k};
        if ~isfield(file_cols, obj), continue; end
        col = file_cols.(obj);
        if ~any(strcmp(col, input_table.Properties.VariableNames)), continue; end
        path = char(string(input_table.(col)(i)));
        if isempty(path) || exist(path, 'file') ~= 2, continue; end
        try
            St = readtable(path, 'TextType', 'string');
        catch
            continue
        end
        tables{end + 1} = St; %#ok<AGROW>
        origins(end + 1) = struct( ...
            'subject', local_get_string(input_table, i, 'subject'), ...
            'run_label', local_get_string(input_table, i, 'run_label'), ...
            'model', row_model, 'object', obj); %#ok<AGROW>
    end
end
end


% =========================================================================
% Utilities
% =========================================================================
function conds = local_filter_conditions(condition_vec, requested)
present = unique(cellstr(string(condition_vec)), 'stable');
if isempty(requested)
    conds = present(:)'; return
end
requested = cellstr(string(requested));
keep = false(size(present));
for i = 1:numel(requested)
    pp = requested{i};
    if contains(pp, '*')
        rx = ['^', regexptranslate('wildcard', pp), '$'];
        hit = ~cellfun('isempty', regexp(present, rx, 'once'));
    else
        hit = strcmp(present, pp);
    end
    keep = keep | hit(:);
end
conds = present(keep)';
end


function s = local_get_string(t, i, col)
s = '';
if any(strcmp(col, t.Properties.VariableNames))
    val = t.(col)(i);
    if isstring(val) || ischar(val), s = char(val);
    elseif iscell(val), s = char(val{1});
    else
        try, s = char(string(val)); catch, end
    end
end
end


function T = local_empty_metrics_table()
T = table('Size', [0 20], ...
    'VariableTypes', {'string','string','string','string','string', ...
        'string','string','string','string','string', ...
        'double','string','double','double','double', ...
        'double','double','double','double','double'}, ...
    'VariableNames', {'subject','run_label','model','object','study_label', ...
        'source','source_kind','source_set','source_name','condition', ...
        'n_lags','reference','peak_lag_seconds','ref_peak_lag_seconds','peak_lag_bias_seconds', ...
        'auc','ref_auc','auc_ratio','shape_corr','shape_r2'});
T.misspec_r2 = double.empty(0, 1);
T.scaled_resid_rmse = double.empty(0, 1);
end


function row = local_metric_row(org, col, skind, sset, sname, cond, reference, m)
row = table( ...
    string(local_org(org, 'subject')), string(local_org(org, 'run_label')), ...
    string(local_org(org, 'model')), string(local_org(org, 'object')), string(""), ...
    string(col), string(skind), string(sset), string(sname), string(cond), ...
    m.n_lags, string(reference), m.peak_lag_seconds, m.ref_peak_lag_seconds, m.peak_lag_bias_seconds, ...
    m.auc, m.ref_auc, m.auc_ratio, m.shape_corr, m.shape_r2, ...
    'VariableNames', {'subject','run_label','model','object','study_label', ...
        'source','source_kind','source_set','source_name','condition', ...
        'n_lags','reference','peak_lag_seconds','ref_peak_lag_seconds','peak_lag_bias_seconds', ...
        'auc','ref_auc','auc_ratio','shape_corr','shape_r2'});
row.misspec_r2 = m.misspec_r2;
row.scaled_resid_rmse = m.scaled_resid_rmse;
end


function s = local_org(org, f)
if isfield(org, f) && ~isempty(org.(f)), s = org.(f); else, s = ''; end
end
