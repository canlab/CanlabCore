function ax = plot_hrf_study_by_subject(study, varargin)
%PLOT_HRF_STUDY_BY_SUBJECT Plot HRFs by subject, averaging duplicate runs.

p = inputParser;
p.addRequired('study', @isstruct);
p.addParameter('Model', 'sfir', @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('SourceModel', '', @(x) ischar(x) || isstring(x));
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x));
p.addParameter('Signature', '', @(x) ischar(x) || isstring(x));
p.addParameter('PlotType', 'fit', @(x) ischar(x) || isstring(x));
p.addParameter('Unit', 'subject', @(x) ischar(x) || isstring(x));
p.addParameter('ShowSE', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('SEAlpha', 0.10, @(x) isscalar(x) && x >= 0 && x <= 1);
p.addParameter('ShowGroupMean', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('ShowGroupSE', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('GroupLineWidth', 2.5, @(x) isscalar(x) && x > 0);
p.addParameter('BaselineSeconds', 0, @(x) isscalar(x) && x >= 0);
p.addParameter('TrialOutlierPolicy', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('TrialOutlierZThreshold', 4, @(x) isscalar(x) && x > 0);
p.addParameter('OutlierWeighting', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('OutlierZThreshold', 4, @(x) isscalar(x) && x > 0);
p.addParameter('CurveWeights', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.parse(study, varargin{:});
opts = p.Results;
opts.Condition = char(opts.Condition);
opts.Signature = char(opts.Signature);
opts.SourceModel = lower(strtrim(char(opts.SourceModel)));
opts.Models = local_model_names(opts.Model);
opts.Model = opts.Models{1};
opts.PlotType = lower(char(opts.PlotType));
opts.Unit = char(opts.Unit);

switch opts.PlotType
    case 'fit'
        if numel(opts.Models) > 1
            ax = local_plot_multiple_models(study, opts);
            return
        end
        [Y_runs, SE_runs, run_subject_ids, run_labels, condition_name, x, skipped] = local_collect_curves(study, opts);
    case {'trialmean', 'trial_mean', 'trials'}
        [Y_runs, SE_runs, run_subject_ids, run_labels, condition_name, x, skipped] = local_collect_trial_means(study, opts);
    otherwise
        error('Unknown PlotType: %s. Use ''fit'' or ''trialmean''.', opts.PlotType);
end
if isempty(Y_runs)
    local_error_no_curves(skipped, study, opts);
end

[Y, SE, subject_ids] = local_aggregate_unit(Y_runs, SE_runs, run_subject_ids, lower(char(opts.Unit)), run_labels);
curve_weights = local_curve_weights(Y, opts);

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
if logical(opts.ShowGroupMean)
    group_mean = local_weighted_mean_omitnan(Y, curve_weights)';
    group_se = local_weighted_sem_omitnan(Y, curve_weights)';
    if logical(opts.ShowGroupSE) && size(Y, 1) > 1 && any(~isnan(group_se))
        fill(ax, [x; flipud(x)], [group_mean + group_se; flipud(group_mean - group_se)], ...
            [0 0 0], 'FaceAlpha', min(opts.SEAlpha * 2, 0.25), ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    plot(ax, x, group_mean, 'k-', 'LineWidth', opts.GroupLineWidth, 'DisplayName', 'Group mean');
end
legend(ax, 'Interpreter', 'none');
title(ax, local_title(opts, condition_name, SE, Y, study, curve_weights), 'Interpreter', 'none');
xlabel(ax, 'Seconds after event onset');
ylabel(ax, local_ylabel(opts, study));
hline(0, 'k-');
end

function local_error_no_curves(skipped, study, opts)
reason_summary = local_skip_reason_summary(skipped);
available = local_available_summary(study);
error(['No valid HRF curves found for Model=%s, SourceModel=%s, Signature=%s, Condition=%s. ' ...
    'Skipped %d result(s). Top reasons: %s. Available: %s'], ...
    char(opts.Model), local_display_or_default(opts.SourceModel), local_display_or_default(opts.Signature), ...
    local_display_or_default(opts.Condition), numel(skipped), reason_summary, available);
end

function txt = local_skip_reason_summary(skipped)
if isempty(skipped)
    txt = 'none recorded';
    return
end
reasons = {skipped.reason};
[u, ~, ic] = unique(reasons, 'stable');
counts = accumarray(ic(:), 1);
n = min(numel(u), 4);
parts = cell(1, n);
for i = 1:n
    parts{i} = sprintf('%s (%d)', u{i}, counts(i));
end
txt = strjoin(parts, '; ');
end

function txt = local_available_summary(study)
parts = {};
if isfield(study, 'source_models') && ~isempty(study.source_models)
    parts{end + 1} = sprintf('source_models={%s}', strjoin(cellstr(string(study.source_models)), ', '));
end
if isfield(study, 'score_names') && ~isempty(study.score_names)
    parts{end + 1} = sprintf('score_names={%s}', strjoin(cellstr(string(study.score_names)), ', '));
end
[models, signatures, conditions] = local_available_from_results(study);
if ~isempty(models)
    parts{end + 1} = sprintf('models={%s}', strjoin(models, ', '));
end
if ~isempty(signatures)
    parts{end + 1} = sprintf('signatures={%s}', strjoin(signatures, ', '));
end
if ~isempty(conditions)
    parts{end + 1} = sprintf('conditions={%s}', strjoin(conditions, ', '));
end
if isempty(parts)
    txt = 'none';
else
    txt = strjoin(parts, '; ');
end
end

function [models, signatures, conditions] = local_available_from_results(study)
models = {};
signatures = {};
conditions = {};
for i = 1:numel(study.results)
    r = study.results{i};
    if isempty(r), continue; end
    if isfield(r, 'fits')
        models = [models; fieldnames(r.fits)]; %#ok<AGROW>
    end
    if isfield(r, 'fits_by_signature') && ~isempty(r.fits_by_signature)
        signatures = [signatures; fieldnames(r.fits_by_signature)]; %#ok<AGROW>
    end
    if isfield(r, 'conditions') && ~isempty(r.conditions)
        conditions = [conditions; cellstr(string(r.conditions(:)))]; %#ok<AGROW>
    end
end
models = unique(models, 'stable');
signatures = unique(signatures, 'stable');
conditions = unique(conditions, 'stable');
end

function txt = local_display_or_default(value)
txt = char(value);
if isempty(txt)
    txt = '<default>';
end
end

function reason = local_result_error(study, idx, fallback)
reason = fallback;
if isfield(study, 'errors') && numel(study.errors) >= idx && ~isempty(study.errors{idx})
    reason = char(string(study.errors{idx}));
end
end

function ax = local_plot_multiple_models(study, opts)
figure; ax = axes; hold(ax, 'on');
colors = lines(numel(opts.Models));
all_skipped = struct('index', {}, 'subject', {}, 'reason', {});
condition_name = char(opts.Condition);
plotted = false(1, numel(opts.Models));

for m = 1:numel(opts.Models)
    model_opts = opts;
    model_opts.Model = opts.Models{m};
    [Y_runs, SE_runs, run_subject_ids, run_labels, this_condition, x, skipped] = local_collect_curves(study, model_opts);
    all_skipped = [all_skipped, skipped]; %#ok<AGROW>
    if isempty(Y_runs)
        continue
    end
    if isempty(condition_name) || strcmp(condition_name, char(opts.Condition))
        condition_name = this_condition;
    end
    [Y, ~, subject_ids] = local_aggregate_unit(Y_runs, SE_runs, run_subject_ids, lower(char(opts.Unit)), run_labels);
    curve_weights = local_curve_weights(Y, model_opts);
    model_color = colors(m, :);
    subject_color = model_color + (1 - model_color) * 0.55;

    if logical(opts.ShowGroupMean)
        group_mean = local_weighted_mean_omitnan(Y, curve_weights)';
        group_se = local_weighted_sem_omitnan(Y, curve_weights)';
        if logical(opts.ShowGroupSE) && size(Y, 1) > 1 && any(~isnan(group_se))
            fill(ax, [x; flipud(x)], [group_mean + group_se; flipud(group_mean - group_se)], ...
                model_color, 'FaceAlpha', min(opts.SEAlpha * 2, 0.20), ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        for s = 1:size(Y, 1)
            plot(ax, x, Y(s, :)', '-', 'Color', subject_color, ...
                'LineWidth', 0.75, 'HandleVisibility', 'off');
        end
        plot(ax, x, group_mean, '-', 'Color', model_color, ...
            'LineWidth', opts.GroupLineWidth, 'DisplayName', opts.Models{m});
    else
        for s = 1:size(Y, 1)
            plot(ax, x, Y(s, :)', '-', 'Color', model_color, ...
                'DisplayName', sprintf('%s | %s', opts.Models{m}, subject_ids{s}));
        end
    end
    plotted(m) = true;
end

if ~any(plotted)
    local_error_no_curves(all_skipped, study, opts);
end

legend(ax, 'Interpreter', 'none');
title(ax, local_multi_model_title(opts, condition_name), 'Interpreter', 'none');
xlabel(ax, 'Seconds after event onset');
ylabel(ax, local_ylabel(opts, study));
hline(0, 'k-');
end

function [Y, SE, subject_ids, run_labels, condition_name, x, skipped] = local_collect_curves(study, opts)
Y = [];
SE = [];
subject_ids = {};
run_labels = {};
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
        skipped = local_skip(skipped, s, subject_id, local_result_error(study, s, 'empty result'), missing_policy);
        continue
    end

    [fit_struct, ok, reason] = local_fit_struct(r, opts.Signature);
    if ~ok
        skipped = local_skip(skipped, s, subject_id, reason, missing_policy);
        continue
    end
    [fit, ok, reason] = local_select_fit(fit_struct, model_name, opts.SourceModel);
    if ~ok
        skipped = local_skip(skipped, s, subject_id, reason, missing_policy);
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
    run_labels{end + 1, 1} = local_run_label(study, s, subject_id); %#ok<AGROW>
end
end

function [Y, SE, subject_ids, run_labels, condition_name, x, skipped] = local_collect_trial_means(study, opts)
Y = [];
SE = [];
subject_ids = {};
run_labels = {};
condition_pattern = char(opts.Condition);
condition_name = condition_pattern;
x = [];
skipped = struct('index', {}, 'subject', {}, 'reason', {});
missing_policy = lower(char(opts.MissingPolicy));

for s = 1:numel(study.results)
    subject_id = local_subject_id(study, s);
    r = study.results{s};
    if isempty(r)
        skipped = local_skip(skipped, s, subject_id, local_result_error(study, s, 'empty result'), missing_policy);
        continue
    end
    if ~isfield(r, 'events') || isempty(r.events)
        skipped = local_skip(skipped, s, subject_id, 'missing events for trialmean plot', missing_policy);
        continue
    end
    if ~isfield(r, 'settings') || ~isfield(r.settings, 'TR') || ~isfield(r.settings, 'window_seconds')
        skipped = local_skip(skipped, s, subject_id, 'missing TR/window_seconds for trialmean plot', missing_policy);
        continue
    end

    [tc, ok, reason] = local_trial_timeseries(r, opts.Signature);
    if ~ok
        skipped = local_skip(skipped, s, subject_id, reason, missing_policy);
        continue
    end

    try
        condition_spec = local_condition_spec(r, condition_pattern);
        if isempty(condition_name) || strcmp(condition_name, condition_pattern)
            condition_name = condition_spec.display_label;
        end
        avg = hrf_average_condition_trials(tc, r.events, condition_spec.matched_conditions, ...
            r.settings.TR, r.settings.window_seconds, ...
            'BaselineSeconds', opts.BaselineSeconds, ...
            'OutlierPolicy', opts.TrialOutlierPolicy, ...
            'OutlierZThreshold', opts.TrialOutlierZThreshold);
    catch err
        skipped = local_skip(skipped, s, subject_id, err.message, missing_policy);
        continue
    end

    y = avg.mean(:)';
    se = avg.sem(:)';
    this_x = avg.time(:);
    if isempty(Y)
        Y = y;
        SE = se;
        x = this_x;
    elseif size(Y, 2) == numel(y)
        if numel(this_x) ~= numel(x) || any(abs(this_x(:) - x(:)) > eps(max(abs(x(:))) + 1))
            skipped = local_skip(skipped, s, subject_id, 'trialmean time axis mismatch', missing_policy);
            continue
        end
        Y(end + 1, :) = y; %#ok<AGROW>
        SE(end + 1, :) = se; %#ok<AGROW>
    else
        skipped = local_skip(skipped, s, subject_id, 'trialmean length mismatch', missing_policy);
        continue
    end
    subject_ids{end + 1, 1} = subject_id; %#ok<AGROW>
    run_labels{end + 1, 1} = local_run_label(study, s, subject_id); %#ok<AGROW>
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

function [Y, SE, subject_ids] = local_aggregate_unit(Y_runs, SE_runs, run_subject_ids, unit, run_labels)
if nargin < 5 || isempty(run_labels)
    run_labels = run_subject_ids;
end
switch unit
    case 'run'
        Y = Y_runs;
        SE = SE_runs;
        subject_ids = run_labels(:);
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

function ttl = local_title(opts, condition_name, SE, Y, study, curve_weights)
sig = char(opts.Signature);
if isempty(sig)
    if local_is_mapscore_study(study)
        sig = 'mean_mapscore';
    else
        sig = 'selected signal';
    end
end
parts = {};
if any(~isnan(SE(:)))
    if any(strcmpi(opts.PlotType, {'trialmean', 'trial_mean', 'trials'}))
        if strcmpi(opts.Unit, 'subject')
            parts{end + 1} = 'within-subject/run trial SEM';
        else
            parts{end + 1} = 'within-run trial SEM';
        end
    else
        if strcmpi(opts.Unit, 'subject')
            parts{end + 1} = 'within-subject/run SE';
        else
            parts{end + 1} = 'within-run SE';
        end
    end
end
if logical(opts.ShowGroupSE) && size(Y, 1) > 1
    parts{end + 1} = 'group SEM';
end
if isempty(parts)
    ribbon_txt = 'ribbon=none (SE unavailable)';
else
    ribbon_txt = sprintf('ribbon=%s', strjoin(parts, ' + '));
end
weight_txt = local_weight_title(opts, curve_weights);
source_model = local_source_model_text(opts);
if any(strcmpi(opts.PlotType, {'trialmean', 'trial_mean', 'trials'}))
    ttl = sprintf('model=trialmean, unit=%s, condition=%s, source=%s, %s%s', ...
        lower(char(opts.Unit)), condition_name, sig, ribbon_txt, weight_txt);
elseif local_is_mapscore_study(study)
    ttl = sprintf('model=%s%s, unit=%s, condition=%s, score=%s, %s%s', ...
        char(opts.Model), source_model, lower(char(opts.Unit)), condition_name, sig, ribbon_txt, weight_txt);
else
    ttl = sprintf('model=%s, unit=%s, condition=%s, source=%s, %s%s', ...
        char(opts.Model), lower(char(opts.Unit)), condition_name, sig, ribbon_txt, weight_txt);
end
end

function ttl = local_multi_model_title(opts, condition_name)
sig = char(opts.Signature);
if isempty(sig)
    sig = 'selected signal';
end
parts = {};
if logical(opts.ShowGroupSE)
    parts{end + 1} = 'group SEM';
end
if logical(opts.ShowGroupMean)
    parts{end + 1} = 'thin lines=subjects/runs';
end
if isempty(parts)
    ribbon_txt = 'ribbon=none';
else
    ribbon_txt = strjoin(parts, ', ');
end
ttl = sprintf('models=%s, unit=%s, condition=%s, source=%s, %s', ...
    strjoin(opts.Models, ' + '), lower(char(opts.Unit)), condition_name, sig, ribbon_txt);
end

function ylab = local_ylabel(opts, study)
if any(strcmpi(opts.PlotType, {'trialmean', 'trial_mean', 'trials'}))
    ylab = 'Observed signal';
    return
end
if strcmpi(opts.Model, 'mapscore') || local_is_mapscore_study(study)
    ylab = 'Pattern expression / map score';
else
    ylab = 'Fitted response amplitude';
end
end

function weights = local_curve_weights(Y, opts)
weights = opts.CurveWeights(:);
if ~isempty(weights)
    if numel(weights) ~= size(Y, 1)
        error('CurveWeights must have one value per plotted %s curve (%d).', ...
            lower(char(opts.Unit)), size(Y, 1));
    end
    weights(~isfinite(weights) | weights < 0) = 0;
    return
end

policy = lower(strtrim(char(opts.OutlierWeighting)));
weights = ones(size(Y, 1), 1);
if strcmp(policy, 'none') || size(Y, 1) < 3
    return
end

z = local_curve_max_abs_robust_z(Y);
threshold = opts.OutlierZThreshold;
switch policy
    case {'exclude', 'omit'}
        weights(z > threshold) = 0;
    case {'huber', 'downweight'}
        weights = min(1, threshold ./ max(z, eps));
    case {'bisquare', 'tukey'}
        u = z ./ threshold;
        weights = (1 - u .^ 2) .^ 2;
        weights(u >= 1) = 0;
    otherwise
        error('Unknown OutlierWeighting: %s. Use none, exclude, huber, or bisquare.', policy);
end
end

function z = local_curve_max_abs_robust_z(Y)
n = size(Y, 1);
center = local_nanmedian(Y, 1);
scale = 1.4826 .* local_nanmedian(abs(Y - repmat(center, n, 1)), 1);
scale_bad = scale == 0 | isnan(scale);
fallback_scale = local_nanstd(Y, 1);
scale(scale_bad) = fallback_scale(scale_bad);
scale(scale == 0 | isnan(scale)) = NaN;
Z = abs((Y - repmat(center, n, 1)) ./ repmat(scale, n, 1));
z = zeros(n, 1);
for i = 1:n
    zi = Z(i, :);
    zi = zi(~isnan(zi));
    if ~isempty(zi)
        z(i) = max(zi);
    end
end
end

function m = local_weighted_mean_omitnan(Y, weights)
m = nan(1, size(Y, 2));
for j = 1:size(Y, 2)
    y = Y(:, j);
    valid = ~isnan(y) & weights > 0;
    if any(valid)
        w = weights(valid);
        m(j) = sum(w .* y(valid)) ./ sum(w);
    end
end
end

function se = local_weighted_sem_omitnan(Y, weights)
se = nan(1, size(Y, 2));
for j = 1:size(Y, 2)
    y = Y(:, j);
    valid = ~isnan(y) & weights > 0;
    if sum(valid) < 2
        continue
    end
    w = weights(valid);
    yy = y(valid);
    mu = sum(w .* yy) ./ sum(w);
    n_eff = (sum(w) .^ 2) ./ sum(w .^ 2);
    var_w = sum(w .* (yy - mu) .^ 2) ./ sum(w);
    se(j) = sqrt(var_w ./ max(n_eff, eps));
end
end

function m = local_nanmedian(X, dim)
if nargin < 2, dim = 1; end
if dim == 1
    m = nan(1, size(X, 2));
    for j = 1:size(X, 2)
        y = X(:, j);
        y = y(~isnan(y));
        if ~isempty(y), m(j) = median(y); end
    end
else
    m = nan(size(X, 1), 1);
    for i = 1:size(X, 1)
        y = X(i, :);
        y = y(~isnan(y));
        if ~isempty(y), m(i) = median(y); end
    end
end
end

function s = local_nanstd(X, dim)
mu = local_nanmean(X, dim);
if dim == 1
    centered = X - repmat(mu, size(X, 1), 1);
    n = sum(~isnan(X), 1);
    centered(isnan(centered)) = 0;
    s = sqrt(sum(centered .^ 2, 1) ./ max(n - 1, 1));
elseif dim == 2
    centered = X - repmat(mu, 1, size(X, 2));
    n = sum(~isnan(X), 2);
    centered(isnan(centered)) = 0;
    s = sqrt(sum(centered .^ 2, 2) ./ max(n - 1, 1));
else
    error('local_nanstd supports dim 1 or 2.');
end
s(n < 2) = NaN;
end

function m = local_nanmean(X, dim)
valid = ~isnan(X);
den = sum(valid, dim);
X(~valid) = 0;
m = sum(X, dim) ./ den;
m(den == 0) = NaN;
end

function txt = local_weight_title(opts, weights)
if strcmpi(char(opts.OutlierWeighting), 'none') && isempty(opts.CurveWeights)
    txt = '';
    return
end
txt = sprintf(', curve_weighting=%s, effective_n=%0.3g', ...
    char(opts.OutlierWeighting), sum(weights));
end

function txt = local_source_model_text(opts)
if isempty(opts.SourceModel)
    txt = '';
else
    txt = sprintf(' source_model=%s', opts.SourceModel);
end
end

function tf = local_is_mapscore_study(study)
tf = (isfield(study, 'source') && contains(char(study.source), 'map_scores')) || ...
    (isfield(study, 'mapscore_success') && any(study.mapscore_success)) || ...
    (isfield(study, 'object') && isfield(study, 'model_name') && strcmpi(char(study.model_name), 'mapscore'));
end

function model_names = local_model_names(model_input)
model_names = cellstr(string(model_input));
model_names = cellfun(@(s) lower(strtrim(s)), model_names, 'UniformOutput', false);
model_names = model_names(~cellfun(@isempty, model_names));
if isempty(model_names)
    error('Model must contain at least one model name.');
end
model_names = unique(model_names, 'stable');
end

function [tc, ok, reason] = local_trial_timeseries(r, signature_name)
tc = [];
ok = false;
reason = '';
if ~isempty(signature_name)
    sig_field = local_signature_struct_field(r, signature_name, 'timeseries_by_signature');
    if ~isempty(sig_field)
        tc = r.timeseries_by_signature.(sig_field);
        ok = true;
        return
    end
    if isfield(r, 'signature_meta') && isfield(r.signature_meta, 'selected_signature') && ...
            strcmp(char(r.signature_meta.selected_signature), signature_name) && ...
            isfield(r, 'timeseries') && ~isempty(r.timeseries)
        tc = r.timeseries;
        ok = true;
        return
    end
    reason = sprintf('missing time series for signature %s', signature_name);
    return
end

if isfield(r, 'timeseries') && ~isempty(r.timeseries)
    tc = r.timeseries;
    ok = true;
else
    reason = 'missing timeseries for trialmean plot';
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
sig_field = local_signature_struct_field(r, sig, 'fits_by_signature');
end

function sig_field = local_signature_struct_field(r, sig, struct_field)
sig_field = '';
if ~isfield(r, struct_field) || isempty(r.(struct_field))
    return
end
S = r.(struct_field);
if isfield(S, sig)
    sig_field = sig;
    return
end
candidate = matlab.lang.makeValidName(sig);
if isfield(S, candidate)
    sig_field = candidate;
    return
end
if isfield(r, 'signature_meta') && isfield(r.signature_meta, 'selected_signatures') && ...
        isfield(r.signature_meta, 'selected_signature_fields')
    names = cellstr(string(r.signature_meta.selected_signatures));
    fields = cellstr(string(r.signature_meta.selected_signature_fields));
    idx = find(strcmp(names, sig), 1);
    if ~isempty(idx) && idx <= numel(fields) && isfield(S, fields{idx})
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

function run_label = local_run_label(study, idx, subject_id)
if isfield(study, 'run_labels') && numel(study.run_labels) >= idx && ~isempty(study.run_labels{idx})
    run_label = char(study.run_labels{idx});
else
    run_label = sprintf('%s_run%02d', subject_id, idx);
end
end
