function ax = plot_hrf_study_by_subject(study, varargin)
%PLOT_HRF_STUDY_BY_SUBJECT Plot HRFs by subject, averaging duplicate runs.

p = inputParser;
p.addRequired('study', @isstruct);
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x));
p.addParameter('Signature', '', @(x) ischar(x) || isstring(x));
p.addParameter('Unit', 'subject', @(x) ischar(x) || isstring(x));
p.addParameter('ShowSE', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('SEAlpha', 0.10, @(x) isscalar(x) && x >= 0 && x <= 1);
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.parse(study, varargin{:});
opts = p.Results;
opts.Condition = char(opts.Condition);
opts.Signature = char(opts.Signature);
opts.Model = char(opts.Model);
opts.Unit = char(opts.Unit);

[Y_runs, SE_runs, run_subject_ids, condition_name, x, skipped] = local_collect_curves(study, opts);
if isempty(Y_runs)
    error('No valid HRF curves found. Skipped %d result(s).', numel(skipped));
end

[Y, SE, subject_ids] = local_aggregate_unit(Y_runs, SE_runs, run_subject_ids, lower(char(opts.Unit)));

figure; ax = axes; hold(ax, 'on');
colors = lines(size(Y, 1));
for s = 1:size(Y, 1)
    y = Y(s, :)';
    se = SE(s, :)';
    if logical(opts.ShowSE) && any(~isnan(se))
        fill(ax, [x; flipud(x)], [y + se; flipud(y - se)], colors(s, :), ...
            'FaceAlpha', opts.SEAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    plot(ax, x, y, 'Color', colors(s, :), 'DisplayName', subject_ids{s});
end
legend(ax, 'Interpreter', 'none');
title(ax, local_title(opts, condition_name, SE), 'Interpreter', 'none');
xlabel(ax, 'Seconds after event onset');
ylabel(ax, local_ylabel(opts));
hline(0, 'k-');
end

function [Y, SE, subject_ids, condition_name, x, skipped] = local_collect_curves(study, opts)
Y = [];
SE = [];
subject_ids = {};
condition_pattern = char(opts.Condition);
condition_name = condition_pattern;
x = [];
skipped = struct('index', {}, 'subject', {}, 'reason', {});
model_name = opts.Model;
missing_policy = lower(char(opts.MissingPolicy));

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
    if ~isfield(fit_struct, model_name)
        skipped = local_skip(skipped, s, subject_id, sprintf('missing model %s', model_name), missing_policy);
        continue
    end

    try
        condition_spec = local_condition_spec(r, condition_pattern);
        if isempty(condition_name) || strcmp(condition_name, condition_pattern)
            condition_name = condition_spec.display_label;
        end
    catch
        skipped = local_skip(skipped, s, subject_id, sprintf('missing condition %s', condition_name), missing_policy);
        continue
    end

    fit = fit_struct.(model_name);
    y = local_mean_omitnan(fit.hrf(:, condition_spec.indices), 2)';
    se = local_fit_se(fit, condition_spec.indices, numel(y));
    this_x = local_fit_time(fit, r, numel(y));
    if isempty(Y)
        Y = y;
        SE = se;
        x = this_x;
    elseif size(Y, 2) == numel(y)
        if numel(this_x) ~= numel(x) || any(abs(this_x(:) - x(:)) > eps(max(abs(x(:))) + 1))
            skipped = local_skip(skipped, s, subject_id, 'HRF time axis mismatch', missing_policy);
            continue
        end
        Y(end + 1, :) = y; %#ok<AGROW>
        SE(end + 1, :) = se; %#ok<AGROW>
    else
        skipped = local_skip(skipped, s, subject_id, 'HRF length mismatch', missing_policy);
        continue
    end
    subject_ids{end + 1, 1} = subject_id; %#ok<AGROW>
end
end

function condition_spec = local_condition_spec(r, condition_name)
if ~isfield(r, 'conditions') || isempty(r.conditions)
    error('missing conditions');
end
specs = hrf_resolve_condition_patterns(r.conditions, condition_name, 'DefaultMode', 'first');
condition_spec = specs(1);
condition_spec = local_apply_condition_group_label(r, condition_spec);
end

function condition_spec = local_apply_condition_group_label(r, condition_spec)
if ~isfield(r, 'condition_groups') || isempty(r.condition_groups) || numel(condition_spec.indices) ~= 1
    return
end
idx = condition_spec.indices(1);
groups = r.condition_groups;
if idx <= numel(groups) && isfield(groups, 'display_label')
    condition_spec.display_label = groups(idx).display_label;
    if isfield(groups, 'matched_conditions')
        condition_spec.matched_conditions = groups(idx).matched_conditions;
    end
end
end

function [Y, SE, subject_ids] = local_aggregate_unit(Y_runs, SE_runs, run_subject_ids, unit)
switch unit
    case 'run'
        Y = Y_runs;
        SE = SE_runs;
        subject_ids = run_subject_ids(:);
    case 'subject'
        subject_ids = unique(run_subject_ids, 'stable');
        Y = nan(numel(subject_ids), size(Y_runs, 2));
        SE = nan(numel(subject_ids), size(Y_runs, 2));
        for i = 1:numel(subject_ids)
            wh = strcmp(run_subject_ids, subject_ids{i});
            Y(i, :) = local_mean_omitnan(Y_runs(wh, :), 1);
            SE(i, :) = local_combine_run_se(Y_runs(wh, :), SE_runs(wh, :));
        end
    otherwise
        error('Unknown Unit: %s. Use ''subject'' or ''run''.', unit);
end
end

function se = local_fit_se(fit, condition_idx, n)
se = nan(1, n);
if isfield(fit, 'se') && ~isempty(fit.se) && all(size(fit.se) == size(fit.hrf))
    se = local_combine_condition_se(fit.hrf(:, condition_idx), fit.se(:, condition_idx))';
end
end

function se = local_combine_condition_se(Y, SE)
n_cond = size(Y, 2);
has_model_se = any(~isnan(SE(:)));
if has_model_se
    valid = ~isnan(SE);
    SE(~valid) = 0;
    n = max(sum(valid, 2), 1);
    se = sqrt(sum(SE .^ 2, 2)) ./ n;
else
    se = nan(size(Y, 1), 1);
end
if n_cond > 1
    condition_sem = local_sem_omitnan(Y, 2);
    if has_model_se
        se = sqrt(se .^ 2 + condition_sem .^ 2);
    else
        se = condition_sem;
    end
end
end

function x = local_fit_time(fit, r, n)
if isfield(fit, 'time') && numel(fit.time) == n
    x = fit.time(:);
elseif isfield(fit, 'lag_seconds') && numel(fit.lag_seconds) == n
    x = fit.lag_seconds(:);
elseif isfield(r, 'settings') && isfield(r.settings, 'TR')
    x = 1 + (0:n - 1)' .* r.settings.TR;
else
    x = (1:n)';
end
end

function se = local_combine_run_se(Y, SE)
n_runs = size(Y, 1);
has_model_se = any(~isnan(SE(:)));
if has_model_se
    valid = ~isnan(SE);
    SE(~valid) = 0;
    n = max(sum(valid, 1), 1);
    se = sqrt(sum(SE .^ 2, 1)) ./ n;
else
    se = nan(1, size(Y, 2));
end
if n_runs > 1
    run_sem = local_sem_omitnan(Y, 1);
    if has_model_se
        se = sqrt(se .^ 2 + run_sem .^ 2);
    else
        se = run_sem;
    end
end
end

function ttl = local_title(opts, condition_name, SE)
sig = char(opts.Signature);
if isempty(sig)
    sig = 'selected signal';
end
ribbon_txt = 'ribbon=none (SE unavailable)';
if any(~isnan(SE(:)))
    if strcmpi(opts.Unit, 'subject')
        ribbon_txt = 'ribbon=within-subject/run SE';
    else
        ribbon_txt = 'ribbon=within-run SE';
    end
end
if strcmpi(opts.Model, 'mapscore')
    ttl = sprintf('model=mapscore, unit=%s, condition=%s, score=%s, %s', ...
        lower(char(opts.Unit)), condition_name, sig, ribbon_txt);
else
    ttl = sprintf('model=%s, unit=%s, condition=%s, source=%s, %s', ...
        char(opts.Model), lower(char(opts.Unit)), condition_name, sig, ribbon_txt);
end
end

function ylab = local_ylabel(opts)
if strcmpi(opts.Model, 'mapscore')
    ylab = 'Pattern expression / map score';
else
    ylab = 'Fitted response amplitude';
end
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
else
    centered = X - repmat(mu, 1, size(X, 2));
end
centered(isnan(centered)) = 0;
n = sum(~isnan(X), dim);
sd = sqrt(sum(centered .^ 2, dim) ./ max(n - 1, 1));
se = sd ./ sqrt(n);
se(n == 0 | n == 1) = NaN;
end

function skipped = local_skip(skipped, idx, subject_id, reason, missing_policy)
if strcmp(missing_policy, 'error')
    error('Skipping is disabled. Result %d (%s): %s', idx, subject_id, reason);
elseif ~strcmp(missing_policy, 'warn') && ~strcmp(missing_policy, 'silent')
    error('Unknown MissingPolicy: %s. Use ''warn'', ''silent'', or ''error''.', missing_policy);
end

skipped(end + 1) = struct('index', idx, 'subject', subject_id, 'reason', reason);
if strcmp(missing_policy, 'warn')
    warning('plot_hrf_study_by_subject:SkippingResult', 'Skipping result %d (%s): %s', idx, subject_id, reason);
end
end

function subject_id = local_subject_id(study, idx)
if isfield(study, 'subject_ids') && numel(study.subject_ids) >= idx
    subject_id = char(study.subject_ids{idx});
else
    subject_id = sprintf('sub-%03d', idx);
end
end
