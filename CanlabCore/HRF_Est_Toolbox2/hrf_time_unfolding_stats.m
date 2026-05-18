function stats = hrf_time_unfolding_stats(study, varargin)
%HRF_TIME_UNFOLDING_STATS Subject-level time-unfolding tests across a study.
% Supports within-subject condition contrasts, signature-specific fits, and
% optional two-group comparisons. Duplicate subject IDs are averaged by
% default before statistical testing.

p = inputParser;
p.addRequired('study', @isstruct);
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('SourceModel', '', @(x) ischar(x) || isstring(x));
p.addParameter('Signature', '', @(x) ischar(x) || isstring(x));
p.addParameter('ConditionA', '', @(x) ischar(x) || isstring(x));
p.addParameter('ConditionB', '', @(x) ischar(x) || isstring(x));
p.addParameter('Group', {}, @(x) iscell(x) || isstring(x));
p.addParameter('Unit', 'subject', @(x) ischar(x) || isstring(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.parse(study, varargin{:});
opts = p.Results;

model_name = char(opts.Model);
source_model = lower(strtrim(char(opts.SourceModel)));
unit = lower(char(opts.Unit));
missing_policy = lower(char(opts.MissingPolicy));
opts.ConditionA = char(opts.ConditionA);
opts.ConditionB = char(opts.ConditionB);
opts.Signature = char(opts.Signature);

[D_runs, run_subject_ids, skipped, matchedA, matchedB] = local_collect_contrasts(study, opts, model_name, source_model, missing_policy);
if isempty(D_runs)
    error('No valid HRF fits found for model "%s". Skipped %d result(s).', model_name, numel(skipped));
end

[D, subject_ids, run_to_subject] = local_aggregate_unit(D_runs, run_subject_ids, unit);
[group_labels, group_labels_by_subject] = local_group_labels(opts.Group, study, run_subject_ids, subject_ids, run_to_subject, unit);

[T, P] = local_ttest_each_timepoint(D);

stats = struct();
stats.time = (1:size(D, 2))';
stats.mean = local_mean_omitnan(D, 1)';
stats.sem = local_sem_omitnan(D, 1)';
stats.p_value = P(:);
stats.t_value = T(:);
stats.significant = P(:) < opts.Alpha;
stats.alpha = opts.Alpha;
stats.conditionA = char(opts.ConditionA);
stats.conditionB = char(opts.ConditionB);
stats.conditionA_matched = matchedA;
stats.conditionB_matched = matchedB;
stats.model = model_name;
stats.source_model = source_model;
stats.signature = char(opts.Signature);
stats.unit = unit;
stats.subject_ids = subject_ids(:);
stats.n_subjects = numel(subject_ids);
stats.n_runs = size(D_runs, 1);
stats.run_subject_ids = run_subject_ids(:);
stats.run_to_subject = run_to_subject(:);
stats.data = D;
stats.run_level_data = D_runs;
stats.skipped = skipped;

if ~isempty(group_labels)
    stats.group = group_labels(:);
    stats.group_by_subject = group_labels_by_subject(:);
    ug = unique(group_labels);
    if numel(ug) == 2
        ix1 = strcmp(group_labels, ug{1});
        ix2 = strcmp(group_labels, ug{2});
        [Tg, Pg] = local_ttest2_each_timepoint(D(ix1, :), D(ix2, :));
        stats.group_labels = ug;
        stats.group_p_value = Pg(:);
        stats.group_t_value = Tg(:);
    end
end
end

function [D, subject_ids, skipped, matchedA, matchedB] = local_collect_contrasts(study, opts, model_name, source_model, missing_policy)
D = [];
subject_ids = {};
skipped = struct('index', {}, 'subject', {}, 'reason', {});
matchedA = {};
matchedB = {};

for s = 1:numel(study.results)
    subject_id = local_subject_id(study, s);
    r = study.results{s};
    if isempty(r)
        skipped = local_skip(skipped, s, subject_id, 'empty result', missing_policy);
        continue
    end

    [fit_struct, ok, reason] = local_fit_struct(r, opts.Signature);
    if ~ok
        skipped = local_skip(skipped, s, subject_id, reason, missing_policy);
        continue
    end
    [fit, ok, reason] = local_select_fit(fit_struct, model_name, source_model);
    if ~ok
        skipped = local_skip(skipped, s, subject_id, reason, missing_policy);
        continue
    end

    h = fit.hrf;
    try
        specA = local_condition_spec(r, opts.ConditionA, 'first');
    catch
        skipped = local_skip(skipped, s, subject_id, sprintf('missing ConditionA %s', opts.ConditionA), missing_policy);
        continue
    end

    y = local_mean_omitnan(h(:, specA.indices), 2)';
    matchedA = local_merge_matches(matchedA, specA.matched_conditions);
    if ~isempty(opts.ConditionB)
        try
            specB = local_condition_spec(r, opts.ConditionB, 'first');
        catch
            skipped = local_skip(skipped, s, subject_id, sprintf('missing ConditionB %s', opts.ConditionB), missing_policy);
            continue
        end
        y = y - local_mean_omitnan(h(:, specB.indices), 2)';
        matchedB = local_merge_matches(matchedB, specB.matched_conditions);
    end

    if isempty(D)
        D = y;
    elseif size(D, 2) == numel(y)
        D(end + 1, :) = y; %#ok<AGROW>
    else
        skipped = local_skip(skipped, s, subject_id, 'HRF length mismatch', missing_policy);
        continue
    end
    subject_ids{end + 1, 1} = subject_id; %#ok<AGROW>
end
end

function [D, subject_ids, run_to_subject] = local_aggregate_unit(D_runs, run_subject_ids, unit)
switch unit
    case 'run'
        D = D_runs;
        subject_ids = run_subject_ids(:);
        run_to_subject = (1:numel(run_subject_ids))';
    case 'subject'
        subject_ids = unique(run_subject_ids, 'stable');
        run_to_subject = zeros(numel(run_subject_ids), 1);
        for r = 1:numel(run_subject_ids)
            run_to_subject(r) = find(strcmp(subject_ids, run_subject_ids{r}), 1);
        end
        D = nan(numel(subject_ids), size(D_runs, 2));
        for i = 1:numel(subject_ids)
            D(i, :) = local_mean_omitnan(D_runs(run_to_subject == i, :), 1);
        end
    otherwise
        error('Unknown Unit: %s. Use ''subject'' or ''run''.', unit);
end
end

function m = local_mean_omitnan(X, dim)
if nargin < 2, dim = 1; end
valid = ~isnan(X);
den = sum(valid, dim);
X(~valid) = 0;
m = sum(X, dim) ./ den;
m(den == 0) = NaN;
end

function se = local_sem_omitnan(X, dim)
if nargin < 2, dim = 1; end
mu = local_mean_omitnan(X, dim);
if dim == 1
    centered = X - repmat(mu, size(X, 1), 1);
elseif dim == 2
    centered = X - repmat(mu, 1, size(X, 2));
else
    error('local_sem_omitnan supports dim 1 or 2.');
end
centered(isnan(centered)) = 0;
n = sum(~isnan(X), dim);
sd = sqrt(sum(centered .^ 2, dim) ./ max(n - 1, 1));
se = sd ./ sqrt(n);
se(n == 0) = NaN;
end

function [group_labels, group_labels_by_subject] = local_group_labels(group_input, study, run_subject_ids, subject_ids, run_to_subject, unit)
group_labels = {};
group_labels_by_subject = {};
if isempty(group_input)
    return
end

g = cellstr(string(group_input));
if numel(g) == numel(study.results)
    valid_run_groups = g(local_valid_run_indices(study, run_subject_ids));
elseif numel(g) == numel(run_subject_ids)
    valid_run_groups = g(:);
elseif numel(g) == numel(subject_ids)
    group_labels = g(:);
    group_labels_by_subject = group_labels;
    return
else
    error('Group labels must match number of study results, valid runs, or tested subjects.');
end

group_labels_by_subject = cell(numel(subject_ids), 1);
for i = 1:numel(subject_ids)
    vals = unique(valid_run_groups(run_to_subject == i), 'stable');
    group_labels_by_subject{i} = vals{1};
end

if strcmp(unit, 'run')
    group_labels = valid_run_groups(:);
else
    group_labels = group_labels_by_subject;
end
end

function ix = local_valid_run_indices(study, run_subject_ids)
ix = false(1, numel(study.results));
remaining = run_subject_ids(:);
for i = 1:numel(study.results)
    sid = local_subject_id(study, i);
    wh = find(strcmp(remaining, sid), 1);
    if ~isempty(wh)
        ix(i) = true;
        remaining(wh) = [];
    end
end
end

function [T, P] = local_ttest_each_timepoint(D)
n_tp = size(D, 2);
P = nan(1, n_tp);
T = nan(1, n_tp);
for t = 1:n_tp
    y = D(:, t);
    y = y(~isnan(y));
    if numel(y) < 2
        continue
    end
    [~, pval, ~, st] = ttest(y);
    P(t) = pval;
    T(t) = st.tstat;
end
end

function [T, P] = local_ttest2_each_timepoint(A, B)
n_tp = size(A, 2);
P = nan(1, n_tp);
T = nan(1, n_tp);
for t = 1:n_tp
    a = A(:, t); b = B(:, t);
    a = a(~isnan(a)); b = b(~isnan(b));
    if numel(a) < 2 || numel(b) < 2
        continue
    end
    [~, pval, ~, st] = ttest2(a, b, 'Vartype', 'unequal');
    P(t) = pval;
    T(t) = st.tstat;
end
end

function spec = local_condition_spec(r, condition_name, default_mode)
if nargin < 3, default_mode = 'first'; end
specs = hrf_resolve_condition_patterns(r.conditions, condition_name, 'DefaultMode', default_mode);
spec = specs(1);
if isfield(r, 'condition_groups') && ~isempty(r.condition_groups) && isscalar(spec.indices)
    idx = spec.indices(1);
    if idx <= numel(r.condition_groups) && isfield(r.condition_groups, 'matched_conditions')
        spec.matched_conditions = r.condition_groups(idx).matched_conditions;
        if isfield(r.condition_groups, 'display_label')
            spec.display_label = r.condition_groups(idx).display_label;
        end
    end
end
end

function out = local_merge_matches(existing, new_matches)
out = unique([existing(:); new_matches(:)], 'stable');
end

function [fit_struct, ok, reason] = local_fit_struct(r, signature_name)
fit_struct = struct();
ok = false;
reason = '';
if isempty(signature_name)
    if isfield(r, 'fits')
        fit_struct = r.fits;
        ok = true;
    else
        reason = 'missing fits';
    end
    return
end

if ~isfield(r, 'fits_by_signature')
    reason = 'missing fits_by_signature';
    return
end

sig_field = local_signature_field(r, signature_name);
if isempty(sig_field)
    reason = sprintf('missing signature %s', signature_name);
    return
end
fit_struct = r.fits_by_signature.(sig_field);
ok = true;
end

function [fit, ok, reason] = local_select_fit(fit_struct, model_name, source_model)
fit = struct();
ok = false;
reason = '';
model_name = lower(char(model_name));
source_model = lower(strtrim(char(source_model)));

if isfield(fit_struct, model_name)
    candidate = fit_struct.(model_name);
    wanted_source = local_requested_source_model(model_name, source_model, candidate);
    if local_fit_matches_source(candidate, wanted_source)
        fit = candidate;
        ok = true;
    else
        reason = sprintf('model %s source model mismatch', model_name);
    end
    return
end

if isfield(fit_struct, 'mapscore')
    candidate = fit_struct.mapscore;
    wanted_source = source_model;
    if isempty(wanted_source) && ~strcmp(model_name, 'mapscore')
        wanted_source = model_name;
    end
    if local_fit_matches_source(candidate, wanted_source)
        fit = candidate;
        ok = true;
        return
    end
end

if isempty(source_model)
    reason = sprintf('missing model %s', model_name);
else
    reason = sprintf('missing model %s for source model %s', model_name, source_model);
end
end

function tf = local_fit_matches_source(fit, source_model)
if isempty(source_model)
    tf = true;
    return
end
if isfield(fit, 'source_model') && ~isempty(fit.source_model)
    tf = strcmpi(char(fit.source_model), source_model);
else
    tf = false;
end
end

function source_model = local_requested_source_model(model_name, source_model, fit)
if isempty(source_model) && local_is_wholebrain_model_name(model_name) && ...
        isfield(fit, 'source_model') && ~isempty(fit.source_model)
    source_model = model_name;
end
end

function tf = local_is_wholebrain_model_name(model_name)
tf = ismember(lower(strtrim(char(model_name))), {'fir', 'sfir', 'canonical', 'spline'});
end

function sig_field = local_signature_field(r, sig)
sig_field = '';
if isfield(r.fits_by_signature, sig)
    sig_field = sig;
    return
end
candidate = matlab.lang.makeValidName(sig);
if isfield(r.fits_by_signature, candidate)
    sig_field = candidate;
    return
end
if isfield(r, 'signature_meta') && isfield(r.signature_meta, 'selected_signatures') && ...
        isfield(r.signature_meta, 'selected_signature_fields')
    names = cellstr(string(r.signature_meta.selected_signatures));
    fields = cellstr(string(r.signature_meta.selected_signature_fields));
    idx = find(strcmp(names, sig), 1);
    if ~isempty(idx) && idx <= numel(fields) && isfield(r.fits_by_signature, fields{idx})
        sig_field = fields{idx};
    end
end
end

function skipped = local_skip(skipped, idx, subject_id, reason, missing_policy)
if strcmp(missing_policy, 'error')
    error('Skipping is disabled. Result %d (%s): %s', idx, subject_id, reason);
elseif ~strcmp(missing_policy, 'warn') && ~strcmp(missing_policy, 'silent')
    error('Unknown MissingPolicy: %s. Use ''warn'', ''silent'', or ''error''.', missing_policy);
end

skipped(end + 1) = struct('index', idx, 'subject', subject_id, 'reason', reason);
if strcmp(missing_policy, 'warn')
    warning('hrf_time_unfolding_stats:SkippingResult', 'Skipping result %d (%s): %s', idx, subject_id, reason);
end
end

function subject_id = local_subject_id(study, idx)
if isfield(study, 'subject_ids') && numel(study.subject_ids) >= idx
    subject_id = char(study.subject_ids{idx});
else
    subject_id = sprintf('sub-%03d', idx);
end
end
